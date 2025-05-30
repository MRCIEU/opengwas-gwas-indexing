import gzip
import logging
import math
import multiprocessing
import os
import pickle
import queue
import threading
import traceback

import redis
import shutil
import subprocess
import time

from collections import defaultdict
from dotenv import load_dotenv
from pysam import VariantFile
from retry import retry

from _oci import OCI


load_dotenv()

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=os.environ['LOGGING_LEVEL'])


oci_instance = OCI()


class GWASIndexing:
    def __init__(self):
        _env = os.environ['ENV']
        self.redis = {
            'tasks': redis.Redis(host=os.environ['REDIS_HOST_' + _env], port=os.environ['REDIS_PORT_' + _env], password=os.environ['REDIS_PASS_' + _env], db=0),
            '5e-8_10000_0.001': redis.Redis(host=os.environ['REDIS_HOST_' + _env], port=os.environ['REDIS_PORT_' + _env], password=os.environ['REDIS_PASS_' + _env], db=1),
            '1e-5_1000_0.8': redis.Redis(host=os.environ['REDIS_HOST_' + _env], port=os.environ['REDIS_PORT_' + _env], password=os.environ['REDIS_PASS_' + _env], db=2),
            'phewas_chr_posalleles': redis.Redis(host=os.environ['REDIS_2_HOST_' + _env], port=os.environ['REDIS_2_PORT_' + _env], password=os.environ['REDIS_2_PASS_' + _env], db=1),
            'phewas_cpallales_gwasid': redis.Redis(host=os.environ['REDIS_2_HOST_' + _env], port=os.environ['REDIS_2_PORT_' + _env], password=os.environ['REDIS_2_PASS_' + _env], db=2)
        }

        self.input_dir_local = os.environ['INPUT_DIR_VCF_' + _env]
        self.input_dir_temp = os.environ['INPUT_DIR_TEMP']
        self.temp_dir = os.environ['TEMP_DIR_BCF']
        self.output_dir = os.environ['OUTPUT_DIR']
        self.output_path_index = f"{self.output_dir}/0_pos_prefix_indices_by_dataset"
        for p in [self.input_dir_temp, self.temp_dir, self.output_dir, self.output_path_index]:
            os.makedirs(p, exist_ok=True)

        self.bcftools = os.environ['BCFTOOLS_BINARY_' + _env]
        self.plink = os.environ['PLINK_BINARY_' + _env]
        self.plink_ref = os.environ['PLINK_REF_' + _env]

        self.chunk_size = 10_000_000
        self.phewas_pval = 0.01

        self.tophits_pval = '5e-8'
        self.tophits_kb = '10000'
        self.tophits_r2 = '0.001'

        self.tophits_wide_pval = '1e-5'
        self.tophits_wide_kb = '1000'
        self.tophits_wide_r2 = '0.8'

    def list_pending_tasks_in_redis(self) -> list:
        """
        Fetch list of GWAS IDs from Redis
        :return: List of GWAS IDs
        """
        tasks = self.redis['tasks'].smembers('tasks_pending')
        logging.info('Number of pending tasks: ' + str(len(tasks)))
        return [t.decode('ascii') for t in tasks]

    def setup(self, gwas_id: str) -> (str, str):
        """
        Create temporary directories for this dataset
        :param gwas_id: GWAS ID
        :return: Temp directory and output directory
        """
        temp_dir = f"{self.temp_dir}/{gwas_id}"
        shutil.rmtree(temp_dir, ignore_errors=True)
        os.makedirs(temp_dir)

        output_dir = f"{self.output_dir}/{gwas_id}"
        shutil.rmtree(output_dir, ignore_errors=True)
        os.makedirs(output_dir)

        return temp_dir, output_dir

    def fetch_files(self, gwas_id: str) -> str:
        """
        Fetch the VCF and TBI files locally, or from OCI
        :param gwas_id: GWAS ID
        :return: Full path to the .vcf.gz file
        """
        vcf_path = f"{self.input_dir_local}/{gwas_id}/{gwas_id}.vcf.gz"
        logging.debug(f"Fetching {gwas_id}")

        if not os.path.exists(vcf_path) or not os.path.exists(vcf_path + '.tbi'):
            dataset_path = f"{self.input_dir_temp}/{gwas_id}"
            shutil.rmtree(dataset_path, ignore_errors=True)
            os.makedirs(dataset_path)
            vcf_path = f"{dataset_path}/{gwas_id}.vcf.gz"
            logging.debug('Downloading', gwas_id)
            with open(vcf_path, 'wb') as f:
                f.write(oci_instance.object_storage_download('data', f"{gwas_id}/{gwas_id}.vcf.gz").data.content)
            with open(f"{vcf_path}.tbi", 'wb') as f:
                f.write(oci_instance.object_storage_download('data', f"{gwas_id}/{gwas_id}.vcf.gz.tbi").data.content)

        return vcf_path

    def get_query_and_print_string(self, vcf_path: str) -> str:
        """
        Get the bcftools query string and awk print string based on availability of SS (sample size)
        :param vcf_path: Full path to the .vcf.gz file
        :return: bcftools_query_string, awk_print_string
        """
        vcf = VariantFile(vcf_path)
        available_columns = next(vcf.fetch()).format.keys()

        if 'SS' in available_columns:
            awk_print_string = f'{{print $1, $2, $3, $4, $5, $6, $7, $8, 10^-$9, $10}}'
        else:
            vcf.seek(0)
            global_fields = [x for x in vcf.header.records if x.key == "SAMPLE"][0]
            if 'TotalControls' in global_fields.keys():
                ss = str(int(float(global_fields['TotalControls'])) + int(float(global_fields.get('TotalCases', 0))))
            else:
                ss = '.'
            awk_print_string = f'{{print $1, $2, $3, $4, $5, $6, $7, $8, 10^-$9, \"{ss}\"}}'

        return '%CHROM %POS %ID %ALT %REF[ %AF %ES %SE %LP %SS]\n', awk_print_string

    def extract_vcf(self, gwas_id: str, vcf_path: str, temp_dir: str, bcftools_query_string: str, awk_print_string: str) -> str:
        """
        Query VCF using bcftools and gzip the output to a temporary file
        :param gwas_id: GWAS ID
        :param vcf_path: Full path to the .vcf.gz file
        :param temp_dir: Temporary directory for this dataset
        :param bcftools_query_string: bcftools query string
        :param awk_print_string: awk print string
        :return: Full path to temporary file produced by bcftools
        """
        query_out_path = f"{temp_dir}/{gwas_id}"
        logging.debug(f"Extracting {gwas_id}")

        cmd = (f"{self.bcftools} query -f'{bcftools_query_string}' {vcf_path}"
               f" | awk '{awk_print_string}'"
               f" | grep -v inf | gzip -c > {query_out_path}")

        if subprocess.call(cmd, shell=True) != 0:
            logging.error(f"FAILED bcftools query {gwas_id}")
            raise Exception(f"FAILED bcftools query {gwas_id}")

        logging.debug(f"Extracted {gwas_id}")
        return query_out_path

    def read_gwas(self, gwas_id: str, query_out_path: str) -> (dict, dict):
        """
        Read the temporary file, store the data in chunks and filter out records that will be used in PheWAS
        :param gwas_id: GWAS ID
        :param query_out_path: Full path to temporary file produced by bcftools
        :return: Dictionary of GWAS data (by chromosome and then position prefix)
            and dictonary of PheWAS data (by chr:pos:ea:nea)
        """
        gwas = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        phewas = {}
        logging.debug(f"Reading bcf {gwas_id}")
        with gzip.open(query_out_path) as f:
            for line in f:
                l = ['' if x == '.' else x for x in line.rstrip().decode('utf-8').split(' ')]
                # gwas[l[0].lstrip('0')][int(l[1]) // self.chunk_size][int(l[1])].append(
                #     [l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9]]
                #     # rsid, ea,  nea,  af,   beta, se,   pval, ss
                # )
                if float(l[8]) < self.phewas_pval:
                    phewas[f"{l[0].lstrip('0')}:{l[1]}:{l[3]}:{l[4]}"] = float(l[8])

        # for chr in gwas.keys():
        #     for pos_prefix in gwas[chr].keys():
        #         gwas[chr][pos_prefix] = dict(sorted(gwas[chr][pos_prefix].items()))
        #     gwas[chr] = dict(sorted(gwas[chr].items()))

        logging.debug(f"Read bcf {gwas_id}")
        return dict(sorted(gwas.items())), phewas

    def write_gwas(self, gwas_id: str, gwas: dict, output_dir: str) -> None:
        """
        Write the GWAS chunks to files and produce the index file
        :param gwas_id: GWAS ID
        :param gwas: Dictionary of GWAS data
        :param output_dir: Output directory for this dataset
        :return: None
        """
        output_path = output_dir
        logging.debug(f"Writing {gwas_id}")

        pos_prefix_index = defaultdict(list)
        for chr in gwas.keys():
            pos_prefix_index[chr] = list(set(gwas[chr].keys()))
            for pos_prefix in gwas[chr].keys():
                with gzip.open(output_path + f"/{chr}_{pos_prefix}", 'wb') as f:
                    pickle.dump(gwas[chr][pos_prefix], f)
        pos_prefix_index = dict(pos_prefix_index)

        with gzip.open(f"{self.output_path_index}/{gwas_id}", 'wb') as f:
            pickle.dump(pos_prefix_index, f)
        logging.debug(f"Written {gwas_id}")

    def save_chunks_and_index(self, gwas_id) -> int:
        """
        Upload the GWAS chunks to OCI and save the index to Redis
        :param gwas_id: GWAS ID
        :return: Number of chunks
        """
        file_list = os.listdir(f"{self.output_dir}/{gwas_id}")
        if not file_list:
            logging.error(f"FAILED No files found for {gwas_id}")
            raise Exception(f"FAILED No files found for {gwas_id}")

        logging.debug(f"Uploading {gwas_id} {len(file_list)}")

        n_workers = int(os.environ['N_PROC']) * 2

        qqueue = queue.Queue()
        for i, filename in enumerate(file_list):
            qqueue.put((i, filename))
        for _ in range(n_workers):
            qqueue.put(None)

        threads = []
        for thread_id in range(n_workers):
            thread = threading.Thread(target=file_upload_worker, args=(gwas_id, thread_id, qqueue))
            thread.start()
            threads.append(thread)
        for thread in threads:
            thread.join()

        oci_instance.object_storage_upload('data-chunks', f"0_pos_prefix_indices_by_dataset/{gwas_id}", open(f"{self.output_path_index}/{gwas_id}", 'rb'))

        with gzip.open(f"{self.output_path_index}/{gwas_id}", 'rb') as f:
            self.redis['tasks'].hset('gwas_pos_prefix_indices', gwas_id, f.read())

        logging.debug(f"Uploaded {gwas_id}")
        return len(file_list)

    def save_phewas(self, gwas_id: str, id_n: str, phewas: dict) -> None:
        """
        Save the PheWAS data to Redis
        :param gwas_id: GWAS ID
        :param id_n: id(n) of the GwasInfo node
        :param phewas: Dictionary of PheWAS data (by chr:pos:ea:nea)
        :return:
        """
        chr_posalleles = defaultdict(list)
        posalleles_pval = {}
        for chr_pos_ea_oa, pval in phewas.items():
            chr, pos, ea, oa = chr_pos_ea_oa.split(':')
            chr_posalleles[chr].append((pos, ea, oa))
            posalleles_pval[chr_pos_ea_oa] = pval
        t = time.time()

        pipe = self.redis['phewas_chr_posalleles'].pipeline()
        for chr, pos_ea_oa_tuples in chr_posalleles.items():
            pipe.zadd(chr, {'{}:{}:{}'.format(peo[0], peo[1], peo[2]): peo[0] for peo in pos_ea_oa_tuples})
        r1 = pipe.execute()
        success_chr_posalleles = sum(1 for r in r1 if r > 0)
        t2 = time.time()

        pipe = self.redis['phewas_cpallales_gwasid'].pipeline()
        for chr_pos_ea_oa, pval in posalleles_pval.items():
            pipe.zadd(chr_pos_ea_oa, {id_n: pval})
        r2 = pipe.execute()
        success_cpallales_gwasid = sum(1 for r in r2 if r > 0)

        logging.debug(f"Saved PheWAS of {gwas_id}: {success_chr_posalleles} ({round(t2 - t, 3)} s), {success_cpallales_gwasid} ({round(time.time() - t2, 3)} s)")

    def extract_vcf_tophits(self, gwas_id: str, vcf_path: str, temp_dir: str, params: dict) -> str:
        """
        Extract tophits from the VCF file using bcftools
        :param gwas_id: GWAS ID
        :param vcf_path: Full path to the .vcf.gz file
        :param temp_dir: Temporary directory for this dataset
        :param params: Parameters used for extracting and clumping tophits
        :return: Full path to temporary (tophits) file produced by bcftools
        """
        tophits_path = f"{temp_dir}/{gwas_id}.tophits"
        logging.debug(f"Extracting tophits {gwas_id}")

        cmd = (f"{self.bcftools} view -i 'FORMAT/LP>{-math.log10(float(params['pval']))}' {vcf_path}"
               f" | {self.bcftools} query -f'%ID[ %ES %SE %LP]\n'"
               f" | awk ' function abs(v) {{return v < 0 ? -v : v}} "
               f'BEGIN {{print "SNP P"}}; {{ if ($2 != 0) print $1, $3/abs($2)}}\''  # header must be p rather than se/beta to be recognised by plink
               f" | awk '!seen[$1]++' > {tophits_path}")

        if subprocess.call(cmd, shell=True) != 0:
            logging.error(f"FAILED bcftools tophits query {gwas_id}")
            raise Exception(f"FAILED bcftools tophits query {gwas_id}")

        logging.debug(f"Extracted tophits {gwas_id} with {params['pval']}_{params['kb']}_{params['r2']}")
        return tophits_path

    def clump_tophits(self, gwas_id: str, tophits_path: str, params: dict) -> str:
        """
        Clump the tophits using plink
        :param gwas_id: GWAS ID
        :param tophits_path: Full path to temporary (tophits) file produced by bcftools
        :param params: Parameters used for extracting and clumping tophits
        :return: Full path to clumped tophits
        """
        suffix = f"{params['pval']}_{params['kb']}_{params['r2']}"
        logging.debug(f"Clumping tophits {gwas_id} for {suffix}")

        cmd = (f"{self.plink} --bfile {self.plink_ref} --clump {tophits_path} "
               f"--clump-kb {params['kb']} --clump-r2 {params['r2']} --clump-p1 1 --clump-p2 1 --out {tophits_path}"
               " > /dev/null")

        if subprocess.call(cmd, shell=True) != 0:
            logging.error(f"FAILED clumping {gwas_id} for {suffix}")
            raise Exception(f"FAILED clumping {gwas_id} for {suffix}")

        clumped_rsids_path = f"{tophits_path}.clumped.{suffix}"
        if os.path.isfile(tophits_path + '.clumped'):
            with open(tophits_path + '.clumped', 'rt') as f, open(clumped_rsids_path, 'wt') as fo:
                next(f)  # Skip the header
                for line in f:
                    line = line.strip("\n").split()
                    if len(line) >= 3:
                        fo.write(f"{line[2]}\n")
        else:
            return ''

        logging.debug(f"Clumped tophits {gwas_id} for {suffix}")
        return clumped_rsids_path

    def extract_vcf_by_rsids(self, gwas_id: str, vcf_path: str, temp_dir: str, rsids_path: str, bcftools_query_string: str, awk_print_string: str, params: dict) -> str:
        """
        Extract tophits from the VCF file using bcftools and rsids
        :param gwas_id: GWAS ID
        :param vcf_path: Full path to the .vcf.gz file
        :param temp_dir: Temporary directory for this dataset
        :param rsids_path: Full path to the file that contains the list of rsids
        :param bcftools_query_string: bcftools query string
        :param awk_print_string: awk print string
        :param params: Parameters used for extracting and clumping tophits
        :return: Full path to temporary (tophits) file produced by bcftools
        """
        suffix = f"{params['pval']}_{params['kb']}_{params['r2']}"
        query_out_path = f"{temp_dir}/{gwas_id}"
        logging.debug(f"Extracting {gwas_id} using rsids for {suffix}")

        cmd = (f"{self.bcftools} view -i 'ID=@{rsids_path}' {vcf_path}"
               f" | {self.bcftools} query -f'{bcftools_query_string}'"
               f" | awk '{awk_print_string}'"
               f" | grep -v inf | gzip -c > {query_out_path}")

        if subprocess.call(cmd, shell=True) != 0:
            logging.error(f"FAILED bcftools query {gwas_id} using rsids")
            raise Exception(f"FAILED bcftools query {gwas_id} using rsids")

        logging.debug(f"Extracted {gwas_id} using rsids for {suffix}")
        return query_out_path

    def read_tophits(self, gwas_id: str, query_out_path: str, params: dict) -> list:
        """
        Read tophits from bcftools query result
        :param gwas_id: GWAS ID
        :param query_out_path: Full path to temporary (tophits) file produced by bcftools
        :param params: Parameters used for extracting and clumping tophits
        :return: List of tophits
        """
        suffix = f"{params['pval']}_{params['kb']}_{params['r2']}"
        gwas = []
        logging.debug(f"Reading bcf {gwas_id} tophits for {suffix}")
        with gzip.open(query_out_path) as f:
            for line in f:
                l = ['' if x == '.' else x for x in line.rstrip().decode('utf-8').split(' ')]
                gwas.append(l[:10])

        gwas = list(sorted(gwas, key=lambda x: float(x[8])))  # Sort by pval, ascending
        logging.debug(f"Read bcf {gwas_id} tophits for {suffix}")
        return gwas

    def write_tophits(self, gwas_id: str, gwas: list, output_dir: str, params: dict) -> None:
        """
        Write the tophits to a temporary file
        :param gwas_id: GWAS ID
        :param gwas: List of tophits
        :param output_dir: Output directory for the tophits file
        :param params: Parameters used for extracting and clumping tophits
        :return:
        """
        suffix = f"{params['pval']}_{params['kb']}_{params['r2']}"
        output_path = f"{output_dir}/tophits.{suffix}"
        logging.debug(f"Writing {gwas_id} tophits {suffix}")

        with gzip.open(output_path, 'wb') as f:
            pickle.dump(gwas, f)
        logging.debug(f"Written {gwas_id} tophits {suffix}")

    def save_tophits(self, gwas_id: str, output_dir: str, params: dict) -> None:
        """
        Upload the tophits to OCI and Redis
        :param gwas_id: GWAS ID
        :param output_dir: Full path to the tophits file
        :param params: Parameters used for extracting and clumping tophits
        :return:
        """
        suffix = f"{params['pval']}_{params['kb']}_{params['r2']}"
        output_path = f"{output_dir}/tophits.{suffix}"
        logging.debug(f"Uploading {gwas_id} tophits for {suffix}")

        with open(output_path, 'rb') as f:
            oci_instance.object_storage_upload('tophits', f"{gwas_id}/{suffix}", f)

        with gzip.open(output_path, 'rb') as f:
            tophits = pickle.loads(f.read())
            members_and_scores = {}
            for t in tophits:
                score = float(t[8])
                del t[8]
                members_and_scores[':'.join(t)] = score
            self.redis[suffix].zadd(gwas_id, members_and_scores)

        logging.debug(f"Uploaded {gwas_id} tophits for {suffix}")

    def cleanup(self, gwas_id: str) -> None:
        """
        Remove temporary files and directories for this dataset
        :param gwas_id: GWAS ID
        :return: None
        """
        logging.debug(f"Cleaning up {gwas_id}")
        shutil.rmtree(f"{self.input_dir_temp}/{gwas_id}", ignore_errors=True)
        shutil.rmtree(f"{self.temp_dir}/{gwas_id}", ignore_errors=True)
        shutil.rmtree(f"{self.output_dir}/{gwas_id}", ignore_errors=True)

    def report_task_status_to_redis(self, gwas_id: str, successful: bool, n_chunks: int) -> None:
        """
        Remove a task from 'tasks_pending' and add it to 'tasks_completed' or 'tasks_failed'
        :param gwas_id: the full GWAS ID
        :param successful: whether the task is successful or not
        :param n_docs: number of Elasticsearch documents
        :return: None
        """
        self.redis['tasks'].srem('tasks_pending', gwas_id)
        if successful:
            self.redis['tasks'].zadd('tasks_completed', {gwas_id: n_chunks})
            logging.info('Reported {} as completed with {} docs'.format(gwas_id, n_chunks))
        else:
            self.redis['tasks'].zadd('tasks_failed', {gwas_id: int(time.time())})
            logging.info('Reported {} as failed'.format(gwas_id))

    def run_for_single_dataset(self, id_n_and_gwas_id: str) -> (bool, int):
        """
        For a single GWAS dataset, fetch files, extract VCF, generate chunk and index files, upload to OCI and cleanup
        :param id_n_and_gwas_id: id(n) and GWAS ID e.g. 315919:ieu-a-2
        :return: task status and number of chunks
        """
        try:
            id_n, gwas_id = id_n_and_gwas_id.split(':', 1)

            vcf_path = self.fetch_files(gwas_id)
            bcftools_query_string, awk_print_string = self.get_query_and_print_string(vcf_path)
            temp_dir, output_dir = self.setup(gwas_id)

            # all associations
            query_out_path = self.extract_vcf(gwas_id, vcf_path, temp_dir, bcftools_query_string, awk_print_string)
            gwas, phewas = self.read_gwas(gwas_id, query_out_path)
            # self.write_gwas(gwas_id, gwas, output_dir)
            # n_chunks = self.save_chunks_and_index(gwas_id)

            # phewas
            n_chunks = -1
            self.save_phewas(gwas_id, id_n, phewas)

            # tophits
            # for params in [{
            #     'pval': self.tophits_pval,
            #     'kb': self.tophits_kb,
            #     'r2': self.tophits_r2
            # }, {
            #     'pval': self.tophits_wide_pval,
            #     'kb': self.tophits_wide_kb,
            #     'r2': self.tophits_wide_r2
            # }]:
            #     tophits_path = self.extract_vcf_tophits(gwas_id, vcf_path, temp_dir, params)
            #     clumped_rsids_path = self.clump_tophits(gwas_id, tophits_path, params)
            #     if clumped_rsids_path != '':  # only proceed if there are significant clump results
            #         query_out_path = self.extract_vcf_by_rsids(gwas_id, vcf_path, temp_dir, clumped_rsids_path, bcftools_query_string, awk_print_string, params)
            #         gwas = self.read_tophits(gwas_id, query_out_path, params)
            #         self.write_tophits(gwas_id, gwas, output_dir, params)
            #         self.save_tophits(gwas_id, output_dir, params)

            self.cleanup(gwas_id)
        except Exception as e:
            logging.error(traceback.format_exc())
            return False, 0
        return True, n_chunks


def gwas_indexing_worker(proc_id: int, queue: multiprocessing.Queue) -> None:
    """
    The worker function for GWAS indexing
    :param proc_id: Process ID
    :param queue: multiprocessing.Queue containing id_n_and_gwas_id (e.g. 315919:ieu-a-2) strings
    :return:
    """
    logging.info(f"Process {proc_id} started")
    gi = GWASIndexing()
    tp = time.time()

    def _run(id_n_and_gwas_id):
        t0 = time.time()
        successful, n_chunks = gi.run_for_single_dataset(id_n_and_gwas_id)
        gi.report_task_status_to_redis(id_n_and_gwas_id, successful, n_chunks)
        logging.info('Task {} completed in {} s'.format(id_n_and_gwas_id, str(round(time.time() - t0, 3))))

    while True:
        task = queue.get()
        if not task:
            logging.info('Process {} ended in {} s'.format(proc_id, str(round(time.time() - tp, 3))))
            break
        i, id_n_and_gwas_id = task
        _run(id_n_and_gwas_id)


def file_upload_worker(gwas_id: str, thread_id: int, queue: queue.Queue) -> None:
    """
    The worker function for uploading chunk files of a single dataset
    :param thread_id: Thread ID
    :param queue: queue.Queue containing gwas_id strings
    :param gwas_id: the current GWAS ID
    :return:
    """
    output_dir = GWASIndexing().output_dir

    @retry(tries=-1, delay=3)
    def _upload(gwas_id, i, filename):
        with open(f"{output_dir}/{gwas_id}/{filename}", 'rb') as f:
            oci_instance.object_storage_upload('data-chunks', f"{gwas_id}/{filename}", f)
        logging.debug(f"Uploading {gwas_id}, chunk {i}")

    while True:
        task = queue.get()
        if not task:
            logging.debug(f"Terminating file upload worker {gwas_id}, {thread_id}")
            break
        i, filename = task
        _upload(gwas_id, i, filename)


if __name__ == '__main__':
    gi = GWASIndexing()

    n_proc = int(os.environ['N_PROC'])
    while True:
        tasks = gi.list_pending_tasks_in_redis()
        # tasks = ['315919:ieu-a-2']
        if len(tasks) > 0:
            mqueue = multiprocessing.Queue()
            for i, id_n_and_gwas_id in enumerate(tasks):
                mqueue.put((i, id_n_and_gwas_id))
            for _ in range(n_proc):
                mqueue.put(None)

            processes = []
            for proc_id in range(n_proc):
                proc = multiprocessing.Process(target=gwas_indexing_worker, args=(proc_id, mqueue))
                proc.start()
                processes.append(proc)
            for proc in processes:
                proc.join()

        time.sleep(60)
