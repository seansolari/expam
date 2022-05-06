#!/usr/bin/env python3
from argparse import ArgumentParser
from multiprocessing import Process
import os
import re
from typing import Iterable, List, Mapping, Set
from ete3 import NCBITaxa
from expam.classify import ResultsPathConfig
from expam.classify.config import load_results_config
from expam.classify.taxonomy import EntrezRequestor, TaxonomyNCBI
from expam.database import FileLocationConfig
from expam.database.config import make_database_config
from expam.utils import ls


def load_expam_taxonomy_data(db_path: str):
    config: FileLocationConfig = make_database_config(db_path)

    tax_obj = TaxonomyNCBI(config)
    taxid_to_lineage = {data[0]: tuple(data[1:]) for data in tax_obj.load_taxid_lineage_map()}
    taxid_to_name = {}
    name_to_taxid = {}

    with open(config.taxon_rank, 'r') as f:
        for line in f:
            try:
                name, taxid, rank = line.strip().split(",")
            except ValueError as e:
                print(line)
                raise e

            taxid_to_name[taxid] = name
            name_to_taxid[name] = taxid

    return taxid_to_lineage, taxid_to_name, name_to_taxid


def find_compatible_names(taxid_to_lineage: Mapping[str, Iterable[str]], taxid_to_name: Mapping[str, str], name_to_taxid: Mapping[str, str], taxids_to_replace: Iterable[str]):
    ete3_tax_obj = NCBITaxa()
    compatible = dict()

    for taxid in taxids_to_replace:
        try:
            name_lineage = taxid_to_lineage[taxid]
        except KeyError:
            name = taxid_to_name[taxid]
            name_lineage = find_lineage(taxid_to_lineage, name)

        for ancestor_name in name_lineage[::-1]:
            ancestor_taxid = name_to_taxid[ancestor_name]

            try:
                ete3_tax_obj.get_lineage(ancestor_taxid)
                compatible[str(taxid)] = str(ancestor_taxid)
                break
            except ValueError:
                continue
        else:
            raise Exception("Can't find compatible taxid for %s." % taxid)

    return compatible


def find_lineage(taxid_to_lineage: dict, taxa_name: str):
    for lineage in taxid_to_lineage.values():
        try:
            index = lineage.index(taxa_name)
            return lineage[:index + 1]
        except ValueError:
            continue
    else:
        raise Exception("Can't find lineage for %s." % taxa_name)


class CleanSummary(Process):
    def __init__(self, file_list: str, bad_taxids_set, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.file_list = file_list
        self.bad_taxids_set = bad_taxids_set

    def run(self):
        for summary_file in self.file_list:
            print("Processing %s..." % summary_file)
            data = []

            with open(summary_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')

                    if parts[0] in self.bad_taxids_set:
                        continue
                    else:
                        data.append(line)

            with open(summary_file, 'w') as f:
                f.write("".join(data))


class CleanRaw(Process):
    def __init__(self, file_list: str, taxid_map, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.file_list = file_list
        self.taxid_map = taxid_map

    def run(self):
        for raw_file in self.file_list:
            print("Processing %s..." % raw_file)
            data = []

            with open(raw_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')

                    try:
                        new_taxid = self.taxid_map[parts[2]]
                        parts[2] = new_taxid
                    except KeyError:
                        pass

                    data.append("\t".join(parts))

            with open(raw_file, 'w') as f:
                f.write("\n".join(data))


def distribute(iterable, n):
    lists = [[] for _ in range(n)]

    for i, item in enumerate(iterable):
        lists[i % n].append(item)

    return lists


def make_summary_results_compatible(bad_taxids: Set[str], output_files: List[str], n: int = 5):
    # Replace names in tax summary files.
    summary_files = distribute(output_files, n)
    procs = [CleanSummary(summary_files[i], bad_taxids) for i in range(n)]

    for proc in procs:
        proc.start()

    for proc in procs:
        proc.join()


def make_raw_results_compatible(taxid_map: Mapping[str, str], output_files: List[str], n: int = 5):    
    # Replace names in tax raw output files.
    raw_files = distribute(output_files, n)
    procs = [CleanRaw(raw_files[i], taxid_map) for i in range(n)]

    for proc in procs:
        proc.start()

    for proc in procs:
        proc.join()


def load_taxids(taxids_file: str):
    with open(taxids_file, 'r') as f:
        taxids = [re.findall(r"\d+", line)[0] for line in f]

    return taxids


def expam_main(db_path: str, results_path: str, bad_taxids_file: str, n: int):
    bad_taxids = load_taxids(bad_taxids_file)
    results_config: ResultsPathConfig = load_results_config(results_path, create=False)
    
    # Map summary files.
    summary_files: List[str] = ls(results_config.tax, ext=".csv")
    make_summary_results_compatible(set(bad_taxids), summary_files, n=n)

    # Map raw output files.
    raw_files: List[str] = ls(results_config.tax_raw, ext=".csv")
    taxid_to_lineage, taxid_to_name, name_to_taxid = load_expam_taxonomy_data(db_path)
    taxid_map = find_compatible_names(taxid_to_lineage, taxid_to_name, name_to_taxid, bad_taxids)
    make_raw_results_compatible(taxid_map, raw_files, n=n)


def load_taxid_cache(cache_path: str) -> Mapping[str, str]:
    with open(cache_path, 'r') as f:
        cache = dict(line.strip().split('\t') for line in f)

    return cache


def cache_taxids(cache_path: str, cache: Mapping[str, str]):
    with open(cache_path, 'w') as f:
        f.write('\n'.join('\t'.join(v) for v in cache.items()))


def retrieve_taxid_information(taxids: Set[str]):
    # Get lineages for each input taxid.
    requestor: EntrezRequestor = EntrezRequestor()
    taxid_lineage, name_taxid_rank = requestor.request_labels("taxonomy", "xml", list(taxids))
    taxid_lineage = {taxid: lineage.split(',') for taxid, lineage in taxid_lineage}

    taxid_to_name = {}
    name_to_taxid = {}
    for name, data in name_taxid_rank.items():
        taxid, rank = data.split(",")

        taxid_to_name[taxid] = name
        name_to_taxid[name] = taxid

    return taxid_lineage, taxid_to_name, name_to_taxid


def kraken_main(results_path: str, bad_taxids_file: str, n: int):
    # Load previously mapped taxids.
    _kraken_taxid_cache = os.path.join(os.path.dirname(bad_taxids_file), 'taxid_cache.txt')
    if os.path.exists(_kraken_taxid_cache):
        taxid_map = load_taxid_cache(_kraken_taxid_cache)
    else:
        taxid_map = {}

    # Get any new taxids.
    bad_taxids = load_taxids(bad_taxids_file)
    new_bad_taxids = {tid for tid in bad_taxids if tid not in taxid_map}
    if new_bad_taxids:
        taxid_lineage, taxid_to_name, name_to_taxid = retrieve_taxid_information(new_bad_taxids)
        new_taxid_map = find_compatible_names(taxid_lineage, taxid_to_name, name_to_taxid, new_bad_taxids)

        taxid_map.update(new_taxid_map)
        cache_taxids(_kraken_taxid_cache, taxid_map)

        # Fix taxids in kraken output.
        kraken_raw_files = ls(results_path, ext='fq')
        assert all(f.endswith("__output.fq") for f in kraken_raw_files)

        # make_raw_results_compatible(taxid_map, kraken_raw_files, n=n)
    else:
        print("No new taxids to fix. Exiting...")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-db', dest='db_path', help="Path to database.")
    parser.add_argument('-o', '--out', dest="results_path", help="Path to results.")
    parser.add_argument('-T', dest="bad_taxids_file", help="Path to file containing taxids to replace.")
    parser.add_argument('-N', default=5, dest="n_procs", help="Number of worker processes to use.")
    parser.add_argument('--kraken', default=False, dest='kraken', action="store_true", help="Process Kraken(2) files.")

    args = parser.parse_args()

    if not args.kraken:
        expam_main(args.db_path, args.results_path, args.bad_taxids_file, int(args.n_procs))
    else:
        kraken_main(args.results_path, args.bad_taxids_file, int(args.n_procs))

    