import collections
import gzip
import logging
import multiprocessing
import os
import pickle
import redis
import shutil
import subprocess
import time

from collections import defaultdict
from dotenv import load_dotenv
from multiprocessing import Process, Queue
from pysam import VariantFile
from retry import retry

from _oci import OCI


load_dotenv()

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=os.environ['LOGGING_LEVEL'])
logging.getLogger('elasticsearch').setLevel(logging.WARNING)


oci_instance = OCI()


class GWASIndexing:
    def __init__(self):
        _env = os.environ['ENV']
        self.redis = redis.Redis(host=os.environ['REDIS_HOST_' + _env], port=os.environ['REDIS_PORT_' + _env], password=os.environ['REDIS_PASS_' + _env], db=os.environ['REDIS_DB_TASKS'])

        self.input_dir_local = os.environ['INPUT_DIR_VCF_' + _env]
        self.input_dir_temp = os.environ['INPUT_DIR_TEMP']
        self.temp_dir = os.environ['TEMP_DIR_BCF']
        self.output_dir = os.environ['OUTPUT_DIR']
        self.output_path_index = f"{self.output_dir}/0_pos_prefix_indices_by_dataset"
        for p in [self.input_dir_temp, self.temp_dir, self.output_dir, self.output_path_index]:
            os.makedirs(p, exist_ok=True)

        self.bcftools = os.environ['BCFTOOLS_BINARY_' + _env]

        self.chunk_size = 10_000_000

    def list_pending_tasks_in_redis(self) -> list:
        """
        Fetch list of GWAS IDs from Redis
        :return: list of GWAS IDs
        """
        self.redis.select(int(os.environ['REDIS_DB_TASKS']))
        tasks = self.redis.smembers('gwas_pending')
        logging.info('Number of pending tasks: ' + str(len(tasks)))
        return [t.decode('ascii') for t in tasks]

    def fetch_files(self, gwas_id: str) -> str:
        """
        Fetch the VCF and TBI files locally, or from OCI
        :param gwas_id: GWAS ID
        :return: Full path to the .vcf.gz file
        """
        vcf_path = f"{self.input_dir_local}/{gwas_id}/{gwas_id}.vcf.gz"
        logging.debug('Fetching', gwas_id)

        if not os.path.exists(vcf_path) or not os.path.exists(vcf_path + '.tbi'):
            vcf_path = f"{self.input_dir_temp}/{gwas_id}"
            shutil.rmtree(vcf_path, ignore_errors=True)
            os.makedirs(vcf_path)
            logging.debug('Downloading', gwas_id)
            with open(f"{vcf_path}/{gwas_id}.vcf.gz", 'wb') as f:
                f.write(oci_instance.object_storage_download('data', f"{gwas_id}/{gwas_id}.vcf.gz").data)
            with open(f"{vcf_path}/{gwas_id}.vcf.gz.tbi", 'wb') as f:
                f.write(oci_instance.object_storage_download('data', f"{gwas_id}/{gwas_id}.vcf.gz.tbi").data)

        return vcf_path

    def extract_vcf(self, gwas_id: str, vcf_path: str) -> str:
        """
        Query VCF using bcftools and gzip the output to a temporary file
        :param gwas_id: GWAS ID
        :param vcf_path: Full path to the .vcf.gz file
        :return: Full path to temporary file produced by bcftools
        """
        temp_path = f"{self.temp_dir}/{gwas_id}"
        shutil.rmtree(temp_path, ignore_errors=True)
        os.makedirs(temp_path)
        temp_path = f"{temp_path}/{gwas_id}"
        logging.debug('Extracting', gwas_id)

        vcf = VariantFile(vcf_path)
        available_columns = next(vcf.fetch()).format.keys()
        if 'SS' in available_columns:
            cmd = (f"{self.bcftools} query -f'%CHROM %POS %ID %ALT %REF[ %AF %ES %SE %LP %SS]\n' {vcf_path}"
                   f" | awk '{{print $1, $2, $3, $4, $5, $6, $7, $8, 10^-$9, $10}}'"
                   f" | grep -v inf | gzip -c > {temp_path}")
        else:
            vcf.seek(0)
            global_fields = [x for x in vcf.header.records if x.key == "SAMPLE"][0]
            if 'TotalControls' in global_fields.keys():
                SS = int(global_fields['TotalControls']) + int(global_fields.get('TotalCases', 0))
            else:
                SS = '.'
            cmd = (f"{self.bcftools} query -f'%CHROM %POS %ID %ALT %REF[ %AF %ES %SE %LP]\n' {vcf_path}"
                   f" | awk '{{print $1, $2, $3, $4, $5, $6, $7, $8, 10^-$9, \"{SS}\"}}'"
                   f" | grep -v inf | gzip -c > {temp_path}")

        ret = subprocess.call(cmd, shell=True)
        if ret != 0:
            logging.error(f"FAILED bcftools query {gwas_id}")
            raise Exception(f"FAILED bcftools query {gwas_id}")
        logging.debug('Extracted', gwas_id)
        return temp_path

    def read_gwas(self, gwas_id: str, temp_path: str) -> dict:
        """
        Read the temporary file and store the data in chunks
        :param temp_path: Full path to temporary file produced by bcftools
        :return: Dictionary of GWAS data, by chromosome and then position prefix
        """
        gwas = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        logging.debug('Reading bcf', gwas_id)
        with gzip.open(temp_path) as f:
            for line in f:
                l = ['' if x == '.' else x for x in line.rstrip().decode('utf-8').split(' ')]
                gwas[l[0].lstrip('0')][int(l[1]) // self.chunk_size][int(l[1])] = [l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9]]

        for chr in gwas.keys():
            for pos_prefix in gwas[chr].keys():
                gwas[chr][pos_prefix] = dict(sorted(gwas[chr][pos_prefix].items()))
            gwas[chr] = dict(sorted(gwas[chr].items()))
        logging.debug('Read bcf', gwas_id)
        return dict(sorted(gwas.items()))

    def write_gwas(self, gwas_id, gwas: dict) -> None:
        """
        Write the GWAS chunks to files and produce the index file
        :param gwas_id: GWAS ID
        :param gwas: Dictionary of GWAS data
        :return: None
        """
        output_path = f"{self.output_dir}/{gwas_id}"
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        logging.debug('Writing', gwas_id)

        pos_prefix_index = collections.defaultdict(set)
        for chr in gwas.keys():
            pos_prefix_index[chr] = set(gwas[chr].keys())
            for pos_prefix in gwas[chr].keys():
                with gzip.open(output_path + f"/{chr}_{pos_prefix}", 'wb') as f:
                    pickle.dump(gwas[chr][pos_prefix], f)

        with gzip.open(f"{self.output_path_index}/{gwas_id}", 'wb') as f:
            pickle.dump(pos_prefix_index, f)
        logging.debug('Written', gwas_id)

    def upload_files(self, gwas_id) -> int:
        """
        Upload the GWAS chunks and the index file to OCI
        :param gwas_id: GWAS ID
        :return: Number of chunks
        """
        file_list = os.listdir(f"{self.output_dir}/{gwas_id}")
        if not file_list:
            logging.error(f"FAILED No files found for {gwas_id}")
            raise Exception(f"FAILED No files found for {gwas_id}")

        logging.debug('Uploading', gwas_id, len(file_list))

        n_workers = int(os.environ['N_PROC']) * 2

        queue = Queue()
        for i, filename in enumerate(file_list):
            queue.put((i, filename))
        for _ in range(n_workers):
            queue.put(None)

        processes = []
        for proc_id in range(n_workers):
            proc = Process(target=file_upload_worker, args=(gwas_id, proc_id, queue))
            proc.start()
            processes.append(proc)
        for proc in processes:
            proc.join()

        oci_instance.object_storage_upload('data-chunks', f"0_pos_prefix_indices_by_dataset/{gwas_id}", open(f"{self.output_path_index}/{gwas_id}", 'rb'))

        logging.debug('Uploaded', gwas_id)
        return len(file_list)

    def cleanup(self, gwas_id: str) -> None:
        """
        Remove temporary files and directories for this dataset
        :param gwas_id: GWAS ID
        :return: None
        """
        logging.debug('Cleaning up', gwas_id)
        shutil.rmtree(f"{self.input_dir_temp}/{gwas_id}", ignore_errors=True)
        shutil.rmtree(f"{self.temp_dir}/{gwas_id}", ignore_errors=True)
        shutil.rmtree(f"{self.output_dir}/{gwas_id}", ignore_errors=True)

    def report_task_status_to_redis(self, gwas_id: str, successful: bool, n_chunks: int) -> None:
        """
        Remove a task from 'gwas_pending' and add it to 'gwas_completed' or 'gwas_failed'
        :param gwas_id: the full GWAS ID
        :param successful: whether the task is successful or not
        :param n_docs: number of Elasticsearch documents
        :return: None
        """
        self.redis.select(int(os.environ['REDIS_DB_TASKS']))
        self.redis.srem('gwas_pending', gwas_id)
        if successful:
            self.redis.zadd('gwas_completed', {gwas_id: n_chunks})
            logging.info('Reported {} as completed with {} docs'.format(gwas_id, n_chunks))
        else:
            self.redis.zadd('gwas_failed', {gwas_id: int(time.time())})
            logging.info('Reported {} as failed'.format(gwas_id))

    def run_for_single_dataset(self, gwas_id) -> (bool, int):
        """
        For a single GWAS dataset, fetch files, extract VCF, generate chunk and index files, upload to OCI and cleanup
        :param gwas_id: GWAS ID
        :return: task status and number of chunks
        """
        try:
            vcf_path = self.fetch_files(gwas_id)
            temp_path = self.extract_vcf(gwas_id, vcf_path)
            gwas = self.read_gwas(gwas_id, temp_path)
            self.write_gwas(gwas_id, gwas)
            n_chunks = self.upload_files(gwas_id)
            self.cleanup(gwas_id)
        except Exception as e:
            logging.debug(e)
            return False, 0
        return True, n_chunks


def gwas_indexing_worker(proc_id: int, queue: multiprocessing.Queue) -> None:
    """
    The worker function for GWAS indexing
    :param proc_id: process ID
    :param queue: multiprocessing.Queue containing gwas_id strings
    :return:
    """
    logging.info(f"Process {proc_id} started")
    gi = GWASIndexing()
    tp = time.time()

    def _run(gwas_id):
        t0 = time.time()
        successful, n_chunks = gi.run_for_single_dataset(gwas_id)
        gi.report_task_status_to_redis(gwas_id, successful, n_chunks)
        logging.info('Task {} completed in {} s'.format(gwas_id, str(round(time.time() - t0, 3))))

    while True:
        task = queue.get()
        if not task:
            logging.info('Process {} ended in {} s'.format(proc_id, str(round(time.time() - tp, 3))))
            break
        i, gwas_id = task
        _run(gwas_id)


def file_upload_worker(gwas_id: str, proc_id: int, queue: multiprocessing.Queue) -> None:
    """
    The worker function for uploading chunk files of a single dataset
    :param proc_id: process ID
    :param queue: multiprocessing.Queue containing gwas_id strings
    :param gwas_id: the current GWAS ID
    :return:
    """
    output_dir = GWASIndexing().output_dir

    @retry(tries=-1, delay=3)
    def _upload(gwas_id, i, filename):
        oci_instance.object_storage_upload('data-chunks', f"{gwas_id}/{filename}", open(f"{output_dir}/{gwas_id}/{filename}", 'rb'))
        logging.info(f"Uploading {gwas_id}, chunk {i}")

    while True:
        task = queue.get()
        if not task:
            logging.debug(f"Terminating file upload worker {gwas_id}, {proc_id}")
            break
        i, filename = task
        _upload(gwas_id, i, filename)


if __name__ == '__main__':
    gi = GWASIndexing()

    n_proc = int(os.environ['N_PROC'])
    while True:
        while len(gwas_ids := gi.list_pending_tasks_in_redis()) > 0:
            queue = Queue()
            for i, gwas_id in enumerate(gwas_ids):
                queue.put((i, gwas_id))
            for _ in range(n_proc):
                queue.put(None)

            processes = []
            for proc_id in range(n_proc):
                proc = Process(target=gwas_indexing_worker, args=(proc_id, queue))
                proc.start()
                processes.append(proc)
            for proc in processes:
                proc.join()

        time.sleep(60)
