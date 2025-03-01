import collections
import gzip
import os
import pickle
import io
import shutil
import time
import uuid

import requests

from dotenv import load_dotenv

from _oci import OCI


load_dotenv()


class GWASQuery:
    def __init__(self):
        self.oci = OCI()

        self.temp_dir = f"{os.environ['TEMP_DIR_QUERY']}/{uuid.uuid4()}"
        os.makedirs(self.temp_dir, exist_ok=True)

        self.chunk_size = 10_000_000

    def __del__(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _convert_sample_size(self, size):
        if size == '':
            return size
        return float(size) if '.' in size else int(size)

    # for each chr, merge the ranges of positions into minimal ranges (https://leetcode.com/problems/merge-intervals/)
    # then work out the groups of chunk to fetch, without knowing which chunks are actually available
    def merge_pos_ranges(self, pos_tuple_list_by_chr: dict[list[tuple[int, int]]]) -> dict[list[int]]:
        merged_pos_by_chr = collections.defaultdict(list)
        for chr in pos_tuple_list_by_chr.keys():
            pos_tuple_list_by_chr[chr].sort(key=lambda x: x[0])
            for pos_tuple in pos_tuple_list_by_chr[chr]:
                if not merged_pos_by_chr[chr] or merged_pos_by_chr[chr][-1][1] < pos_tuple[0]:
                    merged_pos_by_chr[chr].append(pos_tuple)
                else:
                    merged_pos_by_chr[chr][-1] = (merged_pos_by_chr[chr][-1][0], max(merged_pos_by_chr[chr][-1][1], pos_tuple[1]))
            merged_pos_by_chr[chr].sort(key=lambda x: x[0])
        return merged_pos_by_chr

    # filter out the pos prefixes (chunks) available for the given gwas_id, chr and pos_tuple
    def filter_chunks_available(self, pos_prefix_indices: dict, gwas_id: str, chr: str, pos_tuple: tuple[int]) -> set:
        chunks_available = set()
        for pos_prefix in range(pos_tuple[0] // self.chunk_size, pos_tuple[1] // self.chunk_size + 1):
            if chr in pos_prefix_indices[gwas_id] and pos_prefix in pos_prefix_indices[gwas_id][chr]:
                chunks_available.add(pos_prefix)
        return chunks_available

    # fetch the associations available for the given gwas_id, chr and pos prefixes (chunks)
    def fetch_associations_available(self, gwas_id: str, chr: str, pos_prefixes: set) -> dict:
        associations_available = {}
        for pos_prefix in sorted(pos_prefixes):
            chunk_path = f"{self.temp_dir}/{gwas_id}/{chr}_{pos_prefix}"
            local_file_valid = os.path.exists(chunk_path) and os.path.getsize(chunk_path) > 0
            try:
                with gzip.open(chunk_path, 'rb') as f:
                    associations_available = associations_available | pickle.load(f)
            except Exception as e:
                local_file_valid = False
            if not local_file_valid:
                os.makedirs(f"{self.temp_dir}/{gwas_id}", exist_ok=True)
                with open(chunk_path, 'wb+') as f:
                    f.write(self.oci.object_storage_download('data-chunks', f"{gwas_id}/{chr}_{pos_prefix}").data.content)
                with gzip.open(chunk_path, 'rb') as f:
                    associations_available = associations_available | pickle.load(f)
        return associations_available

    # only leave the associations that are within the pos range
    def trim_and_compose_associations(self, gwasinfo: dict, gwas_id: str, chr: str, associations_available: dict, pos_tuple: tuple[int]) -> list:
        pos_available = list(associations_available.keys())
        if not pos_available or pos_tuple[0] > pos_available[-1] or pos_tuple[1] < pos_available[0]:
            return []

        associations = []
        for pos in pos_available:
            if pos_tuple[0] <= pos <= pos_tuple[1]:
                for assoc in associations_available[pos]:
                    associations.append({
                        'id': gwas_id,
                        'trait': gwasinfo[gwas_id]['trait'],
                        'chr': chr,
                        'position': pos,
                        'rsid': assoc[0],
                        'ea': assoc[1],
                        'nea': assoc[2],
                        'eaf': float(assoc[3]) if assoc[3] != '' else '',
                        'beta': float(assoc[4]) if assoc[4] != '' else '',
                        'se': float(assoc[5]) if assoc[5] != '' else '',
                        'p': float(assoc[6]) if assoc[6] != '' else '',
                        'n': self._convert_sample_size(assoc[7])
                    })
        return associations

    def query(self, gwasinfo: dict, gwas_ids: list[str], query: list[str]) -> list:
        pos_tuple_list_by_chr = collections.defaultdict(list)
        for q in query:
            chr, pos = q.split(':')
            pos_start, pos_end = pos.split('-') if '-' in pos else (pos, pos)
            pos_tuple_list_by_chr[chr].append((int(pos_start), int(pos_end)))
        merged_pos_by_chr = self.merge_pos_ranges(pos_tuple_list_by_chr)

        with gzip.GzipFile(fileobj=io.BytesIO(OCI().object_storage_download('data-chunks', '0_pos_prefix_indices').data.content), mode='rb') as f:
            pos_prefix_indices = pickle.loads(f.read())
            print('pos_prefix_indices read')

        result = []
        for gwas_id in gwas_ids:
            for chr in merged_pos_by_chr:
                for pos_tuple in merged_pos_by_chr[chr]:
                    chunks_available = self.filter_chunks_available(pos_prefix_indices, gwas_id, chr, pos_tuple)
                    associations_available = self.fetch_associations_available(gwas_id, chr, chunks_available)
                    associations = self.trim_and_compose_associations(gwasinfo,gwas_id, chr, associations_available, pos_tuple)
                    result.extend(associations)
        return result


if __name__ == '__main__':
    # gq = GWASQuery()
    # gq.collect_pos_prefix_indices()
    # exit()

    gwas_ids = ['ieu-a-2']
    variants = ['7:105561135', '7:105561135-105563135']
    # variants = ['7:105561135', '7:105561135-105563135', '1:1000000-11000500', '1:30000250-13000400', '2:1000000-10000500', '2:21000200-31000800', 'X:32000000-36000000']

    headers = {
        'X-TEST-MODE-KEY': os.environ['BENCHMARK_API_KEY']
    }

    # api
    r0 = requests.post(os.environ['BENCHMARK_API_URL'] + "/associations", data={'id': gwas_ids, 'variant': variants}, headers=headers).json()
    r0 = sorted(r0, key=lambda x: x['position'])

    # pickle
    t0 = time.time()
    gq = GWASQuery()
    gwasinfo = {r['id']: r for r in requests.post(os.environ['BENCHMARK_API_URL'] + "/gwasinfo", data={'id': gwas_ids}, headers=headers).json()}
    r1 = gq.query(gwasinfo, gwas_ids, variants)
    r1 = sorted(r1, key=lambda x: x['position'])
    print(time.time() - t0)

    r0fset = {frozenset(a.items()) for a in r0}
    r1fset = {frozenset(a.items()) for a in r1}

    r0uniq = r0fset - r1fset
    r1uniq = r1fset - r0fset

    anomalies = []

    r0uniq = [{a[0]: a[1] for a in fset} for fset in r0uniq]
    r0uniq_dict = {f"{a['position']}_{a['rsid']}_{a['ea']}_{a['nea']}": a for a in r0uniq}
    if len(r0uniq_dict) != len(r0uniq):
        anomalies.append(['DUPLICATE_ID', 'r0'])

    r1uniq = [{a[0]: a[1] for a in fset} for fset in r1uniq]
    r1uniq_dict = {f"{a['position']}_{a['rsid']}_{a['ea']}_{a['nea']}": a for a in r1uniq}
    if len(r1uniq_dict) != len(r1uniq):
        anomalies.append(['DUPLICATE_ID', 'r1'])

    print(len(r0uniq_dict), len(r1uniq_dict))

    def _compare(dict1, dict2):
        for id in dict1.keys():
            if dict1[id] != dict2[id]:
                diff_keys = {k for k in dict1[id] if dict1[id][k] != dict2[id][k]}
                if diff_keys != {'n'} or int(dict1[id]['n']) != int(dict2[id]['n']):
                    anomalies.append(['DIFF_IN_VALUE', dict1[id], dict2[id], diff_keys])

    _compare(r0uniq_dict, r1uniq_dict)
    _compare(r1uniq_dict, r0uniq_dict)

    assert anomalies == []
