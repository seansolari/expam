#!/usr/bin/env python3
from argparse import ArgumentParser
from multiprocessing import Process
import re
from typing import List
from ete3 import NCBITaxa
from expam.classify import ResultsPathConfig
from expam.classify.config import load_results_config
from expam.classify.taxonomy import TaxonomyNCBI
from expam.database import FileLocationConfig
from expam.database.config import make_database_config
from expam.utils import ls


def find_compatible_names(db_path: str, taxids_to_replace: List[str]):
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


def make_results_compatible(db_path: str, results_config: ResultsPathConfig, bad_taxids: List[str], n: int = 5):
    taxid_map = find_compatible_names(db_path, bad_taxids)
    bad_taxids_set = set(bad_taxids)

    print(taxid_map)

    # Replace names in tax summary files.
    summary_files = distribute(ls(results_config.tax, ext=".csv"), n)
    procs = [CleanSummary(summary_files[i], bad_taxids_set) for i in range(n)]

    for proc in procs:
        proc.start()

    for proc in procs:
        proc.join()
    
    # Replace names in tax raw output files.
    raw_files = distribute(ls(results_config.tax_raw, ext=".csv"), n)
    procs = [CleanRaw(raw_files[i], taxid_map) for i in range(n)]

    for proc in procs:
        proc.start()

    for proc in procs:
        proc.join()


def main(db_path: str, results_path: str, bad_taxids_file: str, n: int):
    with open(bad_taxids_file, 'r') as f:
        bad_taxids = [re.findall(r"\d+", line)[0] for line in f]

    results_config: ResultsPathConfig = load_results_config(results_path, create=False)
    make_results_compatible(db_path, results_config, bad_taxids, n)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-db', dest='db_path', help="Path to database.")
    parser.add_argument('-o', '--out', dest="results_path", help="Path to results.")
    parser.add_argument('-T', dest="bad_taxids_file", help="Path to file containing taxids to replace.")
    parser.add_argument('-N', default=5, dest="n_procs", help="Number of worker processes to use.")

    args = parser.parse_args()
    main(args.db_path, args.results_path, args.bad_taxids_file, int(args.n_procs))

    