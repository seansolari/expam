import datetime
import os
import re
import requests
import shutil
import time
from ctypes import *
from math import floor
from multiprocessing import Array, shared_memory, Value
from os.path import join

import numpy as np
import pandas as pd
from expam.processes import ReadWorker, ClassifyWorker
from expam.run import EXPAM_TEMP_EXT, TIMEOUT, UNION_RATIO, ControlCenter
from expam.sequences import ls
from expam.stores import TablesDb
from expam.tree import Index

#
# Constants.
#

from expam.run import KEYS_DATA_TYPE, \
    KEYS_DATA_TYPE_SIZE, \
    MAP_DATA_TYPE, \
    MAP_DATA_TYPE_SIZE, \
    LCA_MATRIX_PATH, \
    DB_PATH, \
    ACCESSION_ID_PATH

ENTREZ_EMAIL = "ssol0002@student.monash.edu"
ENTREZ_TOOL = "expam"
ENTREZ_API = "8e5edfa1413576bca4f48b4a5520f6295308"
_STEP = 500  # NCBI Database searches at a time.

PHY_RESULTS = "phy"
TAX_RESULTS = "tax"
RAW_RESULTS = "raw"
TEMP_RESULTS = "temp"

CLASSIFICATIONS_NAME = "classifications.csv"
SPLITS_NAME = "splits.csv"
FINAL_NAME = "results.csv"

TAXID_LINEAGE_MAP_NAME = os.path.join("phylogeny", "taxid_lineage.csv")
TAXON_RANK_MAP_NAME = os.path.join("phylogeny", "taxa_rank.csv")

RANK_ORDER = ["0", "root", "top", "superkingdom", "kingdom", "subkingdom", "superclade", "clade", "subclade",
              "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "superorder", "order",
              "suborder", "superfamily", "family", "subfamily", "supergenus", "genus", "subgenus", "species group",
              "species subgroup", "superspecies", "species", "subspecies", "serogroup", "serotype", "superstrain",
              "strain", "substrain", "no rank"]
TAX_RESULTS_HEADER = ["taxid", "%ClassifiedCumulative", "ClassifiedCumulative", "ClassifiedCount",
                      "%SplitCumulative", "SplitCumulative", "SplitCount", "Rank", "Lineage"]

#
# Functions.
#


def prepare_results(n):
    # Create a shared array to put fully classified reads.
    classifications = Array(c_uint32, n, lock=True)
    splits = Array(c_uint32, n, lock=True)
    unclassified = Value(c_uint32, lock=True)

    return classifications, splits, unclassified


def name_to_id(phylogeny_path):
    # Load Index from Newick file.
    with open(phylogeny_path, "r") as f:
        newick_str = f.read().strip()

    _, indexMap = Index.from_newick(newick_str)
    # Define the index dictionary and phi mapping.
    index = {}
    for i, loc in enumerate(indexMap.pool):
        # Define the name mapping.
        index[loc.name] = i

    return index, indexMap


def string_to_tuple(st):
    opener = st.split("(")[1]
    number_string = opener.split(")")[0]
    return tuple(int(i) for i in number_string.split(",") if i != '')


def run_classifier(reads_url, out_dir, k, n, phylogeny_path, keys_shape, values_shape, logging_dir,
                   taxonomy=False, cutoff=0.0, groups=None, keep_zeros=False, cpm=0.0, use_node_names=True,
                   phyla=False, name_taxa=None, colour_list=None, results_name=None, circle_scale=1.0,
                   paired_end=False, alpha=1.0):
    def format_date():
        _date = datetime.datetime.now()
        year, month, day = _date.year, _date.month, _date.day

        date_string = ""
        for part in (year, month, day):
            if len(str(part)) == 1:
                date_string += "0" + str(part)

            else:
                date_string += str(part)

        return date_string

    # Create a dedicated folder to save results.
    # Format of results folders are ./out_dir/results/YYYYMMDD_ID/...
    results_base = os.path.join(out_dir, "results")
    if not os.path.exists(results_base):
        os.mkdir(results_base)

    if results_name is None:
        folder_base = os.path.join(results_base, format_date())

        i = 1
        while True:
            candidate = folder_base + "_" + str(i)

            if not os.path.exists(candidate):
                results_path = candidate
                break

            else:
                i += 1

    else:
        results_path = os.path.join(results_base, results_name)

        if os.path.exists(results_path):
            over_write = input("Results folder %s already exists! Overwrite? (y/n) ")

            if over_write.lower() == "y":
                print("Overwriting directory...")
                shutil.rmtree(results_path)
            else:
                print("Did not overwrite. Exiting...")
                return

    os.mkdir(results_path)
    print("Results directory created at %s." % results_path)

    print("Loading the map and phylogeny.\n")
    # Load the kmer dictionary.
    db_kmers = SharedDB(out_dir, keys_shape, values_shape)

    phy_results_path = os.path.join(results_path, PHY_RESULTS)
    raw_results_path = os.path.join(phy_results_path, RAW_RESULTS)
    temp_results_path = os.path.join(results_path, TEMP_RESULTS)

    for new_path in [phy_results_path, raw_results_path, temp_results_path]:
        os.mkdir(new_path)

    try:
        # Get name mapping.
        index, phylogeny_index = name_to_id(phylogeny_path)

        # Load LCA Matrix.
        lca_matrix = np.load(join(out_dir, LCA_MATRIX_PATH))

        # Initialise the distribution.
        d = Distribution(
            k=k,
            kmer_db=db_kmers,
            index=phylogeny_index,
            lca_matrix=lca_matrix,
            url=reads_url,
            out_dir=results_path,
            logging_dir=logging_dir,
            temp_dir=temp_results_path,
            paired_end=paired_end,
            alpha=alpha
        )

        # Run the classification with n processes.
        d.run(n)

        # Load accession ids.
        with open(join(out_dir, ACCESSION_ID_PATH), 'r') as f:
            sequence_ids = [
                line.strip().split(",")
                for line in f
            ]

        # Update accession ids in map.
        for (fname, accn_id, tax_id) in sequence_ids:
            if accn_id != "None":
                phylogeny_index[fname].accession_id = accn_id

            if tax_id != "None":
                phylogeny_index[fname].taxid = tax_id

        # Output results.
        results = ClassificationResults(
            index=index,
            phylogeny_index=phylogeny_index,
            in_dir=phy_results_path,
            out_dir=results_path,
            groups=groups,
            keep_zeros=keep_zeros,
            cutoff=cutoff,
            cpm=cpm,
            use_node_names=use_node_names,
            phyla=phyla,
            name_taxa=name_taxa,
            colour_list=colour_list,
            circle_scale=circle_scale
        )

        if taxonomy:
            # Attempt to update taxon ids.
            accession_to_taxonomy(out_dir)

            tax_results_path = os.path.join(results_path, TAX_RESULTS)
            os.mkdir(tax_results_path)

            name_to_lineage, taxon_to_rank = load_taxonomy_map(out_dir)
            results.to_taxonomy(name_to_lineage, taxon_to_rank, tax_results_path)

        results.draw_results()

    finally:

        # Close shared memory.
        keys_shm = shared_memory.SharedMemory(
            name=db_kmers.keys_shm_name,
            create=False
        )
        values_shm = shared_memory.SharedMemory(
            name=db_kmers.values_shm_name,
            create=False
        )
        keys_shm.unlink()
        keys_shm.close()
        values_shm.unlink()
        values_shm.close()

        os.rmdir(temp_results_path)


def load_taxonomy_map(out_dir, convert_to_name=True):
    def yield_lines(fh, watch=False):
        for line in fh:
            line = line.strip()
            if line:
                line = line.split(",")
                yield line

    name_tax_id_url = os.path.join(out_dir, ACCESSION_ID_PATH)
    tax_id_lineage_url = os.path.join(out_dir, TAXID_LINEAGE_MAP_NAME)
    taxa_rank_url = os.path.join(out_dir, TAXON_RANK_MAP_NAME)

    # Create map from scientific name --> (taxid, rank).
    taxon_data = {}
    with open(taxa_rank_url, "r") as f:
        for data in yield_lines(f):
            taxon_data[",".join(data[0:-2])] = tuple(data[-2:])

    # Create map from tax_id --> lineage (tuple).
    tax_id_to_lineage = {}
    with open(tax_id_lineage_url, "r") as f:
        for data in yield_lines(f):
            tax_id_to_lineage[data[0]] = tuple(data[1:])

    if not convert_to_name:
        return tax_id_to_lineage, taxon_data

    # Create map from name --> lineage (tuple).
    name_to_lineage = {}
    with open(name_tax_id_url, "r") as f:
        for data in yield_lines(f):
            name_to_lineage[data[0]] = tax_id_to_lineage[data[2]]

    return name_to_lineage, taxon_data


def load_sequence_map(out_dir):
    sequence_data_url = os.path.join(out_dir, ACCESSION_ID_PATH)

    if not os.path.exists(sequence_data_url):
        return []

    with open(sequence_data_url, "r") as f:
        sequences_ids = [
            line.strip().split(",")
            for line in f
        ]

    return sequences_ids


def load_taxid_lineage_map(out_dir):
    map_url = os.path.join(out_dir, TAXID_LINEAGE_MAP_NAME)

    if not os.path.exists(map_url):
        return []

    with open(map_url, "r") as f:
        taxid_map = [
            line.strip().split(",")
            for line in f
        ]

    return taxid_map


def load_rank_map(out_dir):
    rank_url = os.path.join(out_dir, TAXON_RANK_MAP_NAME)

    name_to_rank = {}

    if not os.path.exists(rank_url):
        return name_to_rank

    with open(rank_url, "r") as f:
        for line in f:
            line = line.strip()

            try:
                name_index = line.index(",")
                name_to_rank[line[:name_index]] = line[name_index + 1:]
            except ValueError:
                continue

    return name_to_rank


def to_first_tab(string):
    name = re.findall(r'^[^\t]+(?=\t)', string)

    if not name:
        print(string)
        raise ValueError("No name found!")

    else:
        return name[0]


def accession_to_taxonomy(out_dir):
    """
    Map accession IDs to taxonomic labels.

    :sequence_ids: List - (sequence_id, accession_id, taxon_id)
    :taxon_ranks: List - (taxon_id, taxonomic ranks)
    :taxa_to_rank: Dict - taxon --> rank
    """

    def tuples_to_disk(lst):
        return "\n".join([",".join(item) for item in lst])

    def dict_to_disk(dct):
        return "\n".join([",".join((key, value)) for key, value in dct.items()])

    sequence_ids = load_sequence_map(out_dir)
    taxon_ranks = load_taxid_lineage_map(out_dir)
    taxa_to_rank = load_rank_map(out_dir)

    # Collect taxon ids for unknown organisms.
    accessions_to_be_mapped = []
    taxa_to_be_collected = []

    for (sequence_id, accession_id, taxon_id) in sequence_ids:
        if taxon_id == "None":
            accessions_to_be_mapped.append(accession_id)
        else:
            taxa_to_be_collected.append(taxon_id)

    if accessions_to_be_mapped:
        accession_to_tax = request_tax_ids("nuccore", accessions_to_be_mapped)

        print("Received %d response(s) for ESummary TaxID request!"
              % len(accession_to_tax))

        for i in range(len(sequence_ids)):
            if sequence_ids[i][1] in accession_to_tax:
                sequence_ids[i][2] = accession_to_tax[sequence_ids[i][1]]
                taxa_to_be_collected.append(sequence_ids[i][2])

    # Collect taxonomic lineages for taxa.
    current_taxa = {taxa[0] for taxa in taxon_ranks}
    taxa_to_be_collected = {  # Set so that we collect unique values.
        taxon_id
        for taxon_id in taxa_to_be_collected
        if taxon_id not in current_taxa
    }

    if taxa_to_be_collected:
        taxid_to_taxon, taxon_to_rank = request_labels("taxonomy", "xml", list(taxa_to_be_collected))

        print("Received %d response(s) for EFetch Taxon request!"
              % len(taxid_to_taxon))

        taxon_ranks.extend(taxid_to_taxon)
        taxa_to_rank.update(taxon_to_rank)

    # Save update maps to disk.
    with open(os.path.join(out_dir, ACCESSION_ID_PATH), "w") as f:
        f.write(tuples_to_disk(sequence_ids))

    with open(os.path.join(out_dir, TAXID_LINEAGE_MAP_NAME), "w") as f:
        f.write(tuples_to_disk(taxon_ranks))

    print("Taxonomic lineages written to %s!" % os.path.join(out_dir, TAXID_LINEAGE_MAP_NAME))

    with open(os.path.join(out_dir, TAXON_RANK_MAP_NAME), "w") as f:
        f.write(dict_to_disk(taxa_to_rank))

    print("Taxonomic ranks written to %s!" % os.path.join(out_dir, TAXON_RANK_MAP_NAME))


def request_tax_ids(db, id_list):
    POST_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    taxids, upto = {}, 0

    PARAMS = {
        "tool": ENTREZ_TOOL,
        "email": ENTREZ_EMAIL,
        "api_key": ENTREZ_API,
        "db": db,
    }

    while upto < len(id_list):
        next_requests = id_list[upto:upto + _STEP]

        print("Posting %d UIDs to NCBI Entrez %s."
              % (len(next_requests), db))

        PARAMS["id"] = ",".join(next_requests)
        esummary_request = requests.post(
            url=POST_URL,
            data=PARAMS
        )

        # Parse TaxonIDs from raw results.
        accn_id = tax_id = None

        for line in esummary_request.text.split("\n"):
            if "<DocSum>" in line:  # Start new record.
                accn_id = tax_id = None

            elif 'AccessionVersion' in line:
                accn_id = parse_id(line)

            elif 'TaxId' in line:
                tax_id = parse_tax_id(line)

            elif "</DocSum>" in line:  # Only save complete records.
                if accn_id is not None and tax_id is not None:
                    taxids[accn_id] = tax_id

        upto += _STEP

        time.sleep(1.0)  # Allow server time to breath.

    return taxids


def parse_id(string):
    new_id = re.findall(r'\<Item Name\="AccessionVersion" Type\="String"\>(.*?)\<\/Item\>', string)

    if not new_id:
        raise ValueError("No taxids found!")

    else:
        return new_id[0]


def parse_tax_id(string):
    taxids = re.findall(r'\<Item Name\="TaxId" Type\="Integer"\>(.*?)\<\/Item\>', string)

    if not taxids:
        raise ValueError("No taxids found!")

    else:
        return taxids[0]


def request_labels(db, retmode, id_list):
    POST_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    taxon_labels, ranks, upto = [], {}, 0

    PARAMS = {
        "tool": ENTREZ_TOOL,
        "email": ENTREZ_EMAIL,
        "api_key": ENTREZ_API,
        "db": db,
        "retmode": retmode,
    }

    while upto < len(id_list):
        next_requests = id_list[upto:upto + _STEP]

        print("Posting %d UIDs to NCBI Entrez %s."
              % (len(next_requests), db))

        PARAMS["id"] = ",".join(next_requests)
        efetch_request = requests.post(
            url=POST_URL,
            data=PARAMS
        )

        # Parse taxonomic labels from raw results.
        tax_id = taxa = rank = name = None
        sub_tax_id = sub_rank = sub_name = None
        collect = True

        for line in efetch_request.text.split("\n"):
            if "<Taxon>" in line and collect:  # Start new record.
                tax_id = taxa = name = rank = None

            elif "<Taxon>" in line and not collect:
                sub_tax_id = sub_rank = sub_name = None

            elif "<TaxId>" in line and collect:
                tax_id = parse_single_tax_id(line)

            elif "<TaxId>" in line and not collect:
                sub_tax_id = parse_single_tax_id(line)

            elif "<Lineage>" in line and collect:
                taxa = parse_lineage(line).replace(",", "")

            elif "<Rank>" in line and not collect:
                if sub_name == "cellular organisms":
                    sub_rank = "top"

                else:
                    sub_rank = parse_rank(line)

            elif "<Rank>" in line and collect:
                if name == "cellular organisms":
                    rank = "top"

                else:
                    rank = parse_rank(line)

            elif "<ScientificName>" in line and not collect:
                sub_name = parse_name(line)

            elif "<ScientificName>" in line and collect:
                name = parse_name(line)

            elif "<LineageEx>" in line:
                collect = False

            elif "</Taxon>" in line and not collect:
                ranks[sub_name] = sub_tax_id + "," + sub_rank

            elif "</LineageEx>" in line:
                collect = True

            elif "</Taxon>" in line and collect:
                if tax_id is not None and taxa is not None:
                    lineage = taxa.strip().split("; ")
                    if name not in lineage and name is not None:
                        lineage += [name]

                    taxon_labels.append([tax_id, ",".join(lineage)])

                ranks[name] = tax_id + ',' + rank

        upto += _STEP

        time.sleep(1.0)  # Allow server time to breath.

    return taxon_labels, ranks


def parse_single_tax_id(string):
    taxids = re.findall(r'\<TaxId\>(.*?)\<\/TaxId\>', string)

    if not taxids:
        raise ValueError("No taxids found!")

    else:
        return taxids[0]


def parse_lineage(string):
    lineage = re.findall(r'\<Lineage\>(.*?)\<\/Lineage\>', string)

    if not lineage:
        raise ValueError("No lineages found!")

    else:
        return lineage[0]


def parse_rank(string):
    rank = re.findall(r'\<Rank\>(.*?)\<\/Rank\>', string)

    if not rank:
        raise ValueError("No rank found!")

    else:
        return rank[0]


def parse_name(string):
    name = re.findall(r'\<ScientificName\>(.*?)\<\/ScientificName\>', string)

    if not name:
        raise ValueError("No rank found!")

    else:
        name = re.sub(r"[,]", "", name[0])
        return name


def combine_results(out_dir, taxa_to_rank, cutoff=0):
    final_url = os.path.join(out_dir, FINAL_NAME)

    def where_positive(lst):
        return [x for x, y in enumerate(lst) if y][0]

    def is_number(string):
        return string.lstrip("-").replace(".", "1").replace("e-", "1").isdigit()

    def classified_counter(num):
        return float(num[0])

    def unclassified_counter(lst):
        return len(lst)

    #
    #   Format: [name, rank, classified, unclassified, [...]]
    #

    classified_url = os.path.join(out_dir, TAX_RESULTS, CLASSIFICATIONS_NAME)
    splits_url = os.path.join(out_dir, TAX_RESULTS, SPLITS_NAME)

    methods = {
        "classified": {
            "counter": classified_counter,
            "index": 2,
            "url": classified_url
        },
        "splits": {
            "counter": unclassified_counter,
            "index": 3,
            "url": splits_url
        },
    }

    results = [None, None, 0, 0, []]
    active_ranks = set()

    for method in ("classified", "splits"):
        if not os.path.exists(methods[method]["url"]):
            continue

        with open(methods[method]["url"], "r") as f:
            first_line = f.readline()
            f.seek(0)

            n_cols = first_line.count(",") + 1
            data = [[] for _ in range(n_cols)]

            for line in f:
                row = line.strip().split(",")
                for i in range(len(row)):
                    if row[i] != '':
                        data[i].append(row[i])

            for entry in data:
                check_numbers = [is_number(value) for value in entry]
                count_start = where_positive(check_numbers)

                count = methods[method]["counter"](entry[count_start:])

                if count < cutoff:
                    continue

                level = results

                for taxa in entry[:count_start]:
                    if taxa == "nan":
                        break

                    rank = taxa_to_rank[taxa]

                    if rank not in RANK_ORDER:
                        raise ValueError("Unknown rank %s for taxa %s!" % (rank, taxa))

                    else:
                        active_ranks.add(rank)

                    check_children = [taxa == child[0] for child in level[4]]

                    if not any(check_children):
                        child_index = len(level[4])
                        level[4].append([taxa, rank, 0, 0, []])

                    else:
                        child_index = where_positive(check_children)

                    level[4][child_index][methods[method]["index"]] += count
                    level = level[4][child_index]

    my_rank_order, i = {}, 0
    for rank in RANK_ORDER:

        if rank in active_ranks:
            my_rank_order[rank] = i
            i += 1

    line_format = ["" for _ in range(4 * len(active_ranks))]

    rank_line = line_format.copy()
    for rank in my_rank_order:
        rank_line[4 * my_rank_order[rank]] = rank

    filler = ["", "classified", "unclassified", ""]
    header_line = line_format.copy()
    for i, section in enumerate(line_format):
        header_line[i] = filler[i % len(filler)]

    formatted_results = ",".join(rank_line) + "\n" + ",".join(header_line) + "\n"

    def format_lineage(level, line):
        my_start_ind = my_rank_order[level[1]] * 4

        line[my_start_ind] = level[0]
        line[my_start_ind + 1] = str(level[2])
        line[my_start_ind + 2] = str(level[3])

        if len(level[4]) > 0:
            lines = [line.copy(), ] + [line_format.copy() for _ in range(len(level[4]) - 1)]

            return "\n".join([format_lineage(level[4][c], lines[c]) for c in range(len(lines))])

        else:
            return ",".join(line)

    formatted_results += format_lineage(results[4][0], line_format.copy())

    with open(final_url, "w") as f:
        f.write(formatted_results)

    print("Formatted results written to %s!" % final_url)


#
# Classes.
#


class SharedDB:
    def __init__(self, out_dir, keys_shape, values_shape):
        #
        print("Preparing shared memory allocations...")

        # Prepare shm segment for keys.
        keys_shm = shared_memory.SharedMemory(
            size=np.product(keys_shape) * KEYS_DATA_TYPE_SIZE,
            create=True
        )
        keys_interface = np.ndarray(
            shape=keys_shape,
            dtype=KEYS_DATA_TYPE,
            buffer=keys_shm.buf
        )

        # Prepare shm segment for values.
        values_shm = shared_memory.SharedMemory(
            size=np.product(values_shape) * MAP_DATA_TYPE_SIZE,
            create=True
        )
        values_interface = np.ndarray(
            shape=values_shape,
            dtype=MAP_DATA_TYPE,
            buffer=values_shm.buf
        )

        db = TablesDb(join(out_dir, DB_PATH))

        #
        print("Loading database keys...")
        db.read_keys(arr=keys_interface)

        #
        print("Loading database values...")
        db.read_values(arr=values_interface)

        db.close()

        # Remember the names for attaching in future.
        self.keys_shm_name = keys_shm.name
        self.values_shm_name = values_shm.name

        # Remember array data.
        self.keys_meta = {
            "shape": keys_shape,
            "dtype": KEYS_DATA_TYPE
        }
        self.values_meta = {
            "shape": values_shape,
            "dtype": MAP_DATA_TYPE
        }

        # Close my access.
        keys_shm.close()
        values_shm.close()


class ExpamClassifierProcesses(ControlCenter):
    @classmethod
    def from_method_dict(cls, logging_dir, config):
        base_center = super(ExpamClassifierProcesses, cls).from_dict(logging_dir, config)

        base_arguments = {
            "group_name": "group_name",
            "workers": "workers",
            "child_statuses": "_child_statuses",
            "phases": "phases",
            "phase_queues": "phase_queues",
            "child_queues": "child_queues",
            "children": "_children",
            "processors": "_processors",
            "transitions": "_transitions",
            "timeout": "_timeout",
        }
        base_attributes = {
            attr: getattr(base_center, attr_reference)
            for attr, attr_reference in base_arguments.items()
        }

        # Create child class from this instance.
        control_center = ExpamClassifierProcesses(logging_dir=logging_dir, **base_attributes)

        control_center.set_methods(
            {
                "classify": {
                    "processor": None,
                    "transition": None
                }
            }
        )

        return control_center


class Distribution:
    def __init__(self, k, kmer_db, index, lca_matrix, url, out_dir, temp_dir, logging_dir, alpha,
                 keep_zeros=False, cutoff=0.0, cpm=0.0, paired_end=False):
        """

        :param k: Int
        :param kmer_db: SharedDB
        :param lca_matrix: np.ndarray
        :param url: str
        :param out_dir: str
        """
        # Set up processing references.
        self.k = k
        self.alpha = alpha
        self.kmer_db = kmer_db
        self.lca_matrix = lca_matrix
        self.out_dir = out_dir
        self.temp_dir = temp_dir
        self.logging_dir = logging_dir

        self.index = index
        self.node_names = [node.name if i > 0 else "unclassified" for i, node in enumerate(index.pool)]

        self.keep_zeros = keep_zeros
        self.cutoff = cutoff
        self.cpm = cpm

        # Url containing reads.
        self.url = url
        self.paired_end = paired_end

    def run(self, n):
        """
        Begin processing reads in self.url.
        We have some collection of processes accessing reads from storage
        and then another collection of processes that analyse the reads.

        """
        file_queue = self.prepare_queue(self.paired_end)
        n_classifiers = floor(UNION_RATIO * n)
        n_readers = n - n_classifiers

        mp_config = {
            "name": "expam",
            "phases": ["classify"],
            "layers": {
                "class": ReadWorker,
                "class_args": {
                    "paired_mode": self.paired_end,
                    "n_classifiers": n_classifiers,
                    "logging_dir": self.logging_dir,
                },
                "parent": ControlCenter,
                "child": {
                    "class": ClassifyWorker,
                    "class_args": {
                        "k": self.k,
                        "kmer_db": self.kmer_db,
                        "lca_matrix": self.lca_matrix,
                        "out_dir": join(self.out_dir, PHY_RESULTS),
                        "temp_dir": self.temp_dir,
                        "n_readers": n_readers,
                        "logging_dir": self.logging_dir,
                        "alpha": self.alpha
                    },
                    "parent": ReadWorker,
                    "child": None,
                    "n": n_classifiers
                },
                "n": n_readers
            },
            "timeout": TIMEOUT
        }

        process_network = ExpamClassifierProcesses.from_method_dict(self.logging_dir, mp_config)

        # Load jobs.
        for file_dir in file_queue:
            process_network.load_job(phase="classify", job=file_dir)

        process_network.run()

        # Combine raw read output.
        temporary_files = ls(self.temp_dir, ext=".reads_%s" % EXPAM_TEMP_EXT)
        base_file_names = {
            re.match(r'(\S+)_\d+.reads_%s' % EXPAM_TEMP_EXT, file_name).group(1)
            for file_name in temporary_files
        }

        result_files = []
        for base_file in base_file_names:
            results_file_name = os.path.join(self.out_dir, PHY_RESULTS, RAW_RESULTS, base_file + ".csv")
            results = [
                self.get_data(join(self.temp_dir, file))
                for file in temporary_files
                if file[:len(base_file)] == base_file
            ]

            with open(results_file_name, "w") as f:
                f.write("\n".join(results))

            result_files.append((results_file_name, base_file))

        # Count reads at each node, for each sample.
        counts_dicts = {
            "C": {
                "out_dir": os.path.join(self.out_dir, PHY_RESULTS, "classified_counts.csv"),
                "counts": {}
            },
            "S": {
                "out_dir": os.path.join(self.out_dir, PHY_RESULTS, "splits_counts.csv"),
                "counts": {}
            }
        }
        counts_dicts["U"] = counts_dicts["C"]

        for result_file, sample_name in result_files:
            counts_dicts["C"]["counts"][sample_name] = [0 for _ in range(len(self.node_names))]
            counts_dicts["S"]["counts"][sample_name] = [0 for _ in range(len(self.node_names))]

            with open(result_file, "r") as f:
                for line in f:
                    data = line.split("\t")
                    counts_dicts[data[0]]["counts"][sample_name][int(data[2].lstrip("p"))] += 1

        # Get non-zero rows.
        non_zero_indices = list({
            i
            for _, sample_name in result_files
            for result_type in ("C", "S")
            for i, count in enumerate(counts_dicts[result_type]["counts"][sample_name])
            if (count > 0) or (i == 0)
        })

        def get_cumulative(lst):
            cumulative_counts = lst.copy()

            for i in range(len(self.node_names) - 1, 0, -1):
                cumulative_counts[i] += sum([
                    cumulative_counts[self.index._pointers[child]]
                    for child in self.index.pool[i].children
                ])

            return cumulative_counts

        # Sample-wise summary of counts in kraken-like format.
        summary_format = "%f%%\t%d\t%d"

        # Format phylogenetic node string.
        self.node_names = [
            "p" + self.node_names[i]
            if self.index.pool[i].type == "Branch" and i > 0
            else self.node_names[i]
            for i in range(len(self.node_names))
        ]

        for _, sample_name in result_files:
            total_reads = sum(counts_dicts["C"]["counts"][sample_name]) + sum(counts_dicts["S"]["counts"][sample_name])

            cumulative_classifications = get_cumulative(counts_dicts["C"]["counts"][sample_name])
            cumulative_splits = get_cumulative(counts_dicts["S"]["counts"][sample_name])

            # Employ cutoff.
            if self.cutoff > 0 or self.cpm > 0:
                cutoff = max((total_reads / 1e6) * self.cpm, self.cutoff)
                valid_inds = list({
                    i
                    for cumulative_count in (cumulative_classifications, cumulative_splits)
                    for i, count in enumerate(cumulative_count)
                    if (count >= cutoff) or (i == 0)
                })
                valid_inds.sort()

            else:
                valid_inds = list(range(len(self.node_names)))

            sample_summary = ""
            for i in valid_inds:
                if i not in non_zero_indices and not self.keep_zeros:
                    continue

                if i == 0:
                    line_summary = "unclassified\t" + summary_format % (
                        round(100 * cumulative_classifications[0] / total_reads, 3),
                        cumulative_classifications[0],
                        counts_dicts["C"]["counts"][sample_name][0]) + ("\t" * 3) + "\n"
                else:
                    line_classified_summary = self.node_names[i] + "\t" + summary_format % (
                        round(100 * cumulative_classifications[i] / total_reads, 3),
                        cumulative_classifications[i],
                        counts_dicts["C"]["counts"][sample_name][i])
                    line_splits_summary = summary_format % (
                        round(100 * cumulative_splits[i] / total_reads, 3),
                        cumulative_splits[i],
                        counts_dicts["S"]["counts"][sample_name][i])
                    line_summary = line_classified_summary + "\t" + line_splits_summary + "\n"

                sample_summary += line_summary

            with open(os.path.join(self.out_dir, PHY_RESULTS, sample_name + ".csv"), "w") as f:
                f.write(sample_summary)

        # Save raw counts to disk.
        do_unclassified = True
        for counts_dict in (counts_dicts["C"], counts_dicts["S"]):
            out_dir, counts = counts_dict["out_dir"], counts_dict["counts"]
            df = pd.DataFrame.from_dict(counts)

            if not self.keep_zeros:
                df = df.iloc[non_zero_indices]
                index = [self.node_names[i] for i in non_zero_indices]
            else:
                index = self.node_names

            if "unclassified" in index and not do_unclassified:
                df.drop(0, inplace=True)
                index = index[1:]

            df.index = index
            df.to_csv(out_dir, sep="\t")

            do_unclassified = False

    @staticmethod
    def get_data(file_url):
        with open(file_url, "r") as f:
            data = f.read().strip()

        # Delete temporary file.
        os.remove(file_url)

        return data

    @staticmethod
    def best_match(file_name, file_list):
        def string_intersect(str_one, str_two):
            l = min(len(str_one), len(str_two))

            for i in range(l):
                if str_one[i] != str_two[i]:
                    return str_one[:i]
            else:
                return str_one[:l]

        match_ind = 0
        best_match = string_intersect(file_name, file_list[0])

        for i in range(1, len(file_list)):
            current_match = string_intersect(file_name, file_list[i])

            if len(current_match) > len(best_match):
                best_match = current_match
                match_ind = i

        return match_ind

    def prepare_queue(self, paired_end=False):
        # Create a list of all the files containing reads to be processed.
        file_list = [
            file_name
            for file_name in ls(self.url)
            if re.match(r"(\S+)\.(?:fa|fna|ffn|fq|fasta|fastq)(?:\.tar)*(?:\.gz)*$", file_name)
        ]

        file_queue = []

        if paired_end:
            while file_list:
                next_file = file_list.pop()

                # Find buddy.
                match_ind = self.best_match(next_file, file_list)
                paired_file = file_list.pop(match_ind)

                file_queue.append((join(self.url, next_file), join(self.url, paired_file)))

        else:
            for file_name in file_list:
                file_queue.append((join(self.url, file_name),))

        return file_queue

    @staticmethod
    def _initialise_results(names, unclassified=True):
        """
        Results are returned in spectral form:
            coord --> # hits at coord.

        :param names: List - node names in the tree.
        :return: pd.Series()

        """
        index = (unclassified * ["unclassified"]) + names
        return pd.Series(
            data=np.zeros(len(names) + unclassified, dtype=int),
            index=index)


class ClassificationResults:
    def __init__(self, index, phylogeny_index, in_dir, out_dir, groups=None, keep_zeros=False, cutoff=0.0, cpm=0.0,
                 use_node_names=True, phyla=False, name_taxa=None, colour_list=None, circle_scale=1.0):
        self.index = index
        self.phylogeny_index = phylogeny_index

        self.in_dir = in_dir
        self.out_dir = out_dir

        self.groups = groups  # [(colour, (name1, name2, ...)), (colour, (name3, ...)), ...]
        self.keep_zeros = keep_zeros

        self.cutoff = cutoff  # Cutoff by counts.
        self.cpm = cpm  # Cutoff by counts per million.
        self.use_node_names = use_node_names

        self.phyla = phyla
        self.colour_list = colour_list
        self.name_taxa = name_taxa
        self.circle_scale = circle_scale

        if self.phyla and self.name_taxa is None:
            raise Exception("Inclusion of phyla requires mapping of names to taxonomic lineages!")

        self.tax_id_hierarchy = {"1": set()}  # Map from tax_id -> immediate children.
        self.tax_id_pool = ["1"]  # Children must appear later than parent this list.

    def to_taxonomy(self, name_to_lineage, taxon_to_rank, tax_dir):
        col_names = ["c_perc", "c_cumul", "c_count", "s_perc", "s_cumul", "s_count", "rank", "scientific name"]

        raw_counts_dir = os.path.join(self.in_dir, 'raw')
        raw_output_dir = os.path.join(tax_dir, 'raw')

        if not os.path.exists(raw_output_dir):
            os.mkdir(raw_output_dir)

        class_counts = pd.read_csv(os.path.join(self.in_dir, "classified_counts.csv"), sep="\t", index_col=0, header=0)
        split_counts = pd.read_csv(os.path.join(self.in_dir, "splits_counts.csv"), sep="\t", index_col=0, header=0)

        # Get rid of phylogenetically printed node names.
        def fix_index(df):
            df.index = [
                index.lstrip("p")
                if index not in self.phylogeny_index._pointers
                else index
                for index in df.index
            ]

        fix_index(class_counts)
        fix_index(split_counts)

        samples = class_counts.columns.tolist()
        phylogeny_nodes = split_counts.index.tolist()

        # Map phylogeny index to taxonomy.
        phylogeny_to_taxonomy, tax_id_data = {}, {}
        for phylogeny_node in phylogeny_nodes:
            tax_id, tax_id_map = self.map_phylogeny_node(phylogeny_node, name_to_lineage, taxon_to_rank)

            phylogeny_to_taxonomy[phylogeny_node] = tax_id
            tax_id_data.update(tax_id_map)

        # Get relevant taxon ids.
        tax_ids, ranks, lineages = ["1"], ["root"], [""]
        for tax_id in self.tax_id_pool[1:]:
            rank, lineage = tax_id_data[tax_id]

            tax_ids.append(tax_id)
            ranks.append(rank)
            lineages.append(lineage)

        for sample_name in samples:
            #
            # Map summary output to taxonomy.
            #

            df = pd.DataFrame(0, index=["unclassified"] + tax_ids, columns=col_names)
            total_counts = sum(class_counts[sample_name]) + sum(split_counts[sample_name])

            df.loc[tax_ids[0]:, "rank"] = ranks
            df.loc[tax_ids[0]:, "scientific name"] = lineages

            # Manually insert unclassified counts.
            df.loc["unclassified", "c_count"] = class_counts.loc["unclassified", sample_name]
            df.loc["unclassified", "c_cumul"] = df.loc["unclassified", "c_count"]
            df.loc["unclassified", "c_perc"] = round(100 * df.loc["unclassified", "c_cumul"] / total_counts, 3)

            for col_prefix, counts_df in (("c_", class_counts), ("s_", split_counts)):
                for node in counts_df.index[1:]:
                    df.loc[phylogeny_to_taxonomy[node], col_prefix + "count"] += counts_df.loc[node, sample_name]

            # Get cumulative counts.
            for tax_id in self.tax_id_pool[::-1]:
                for pre in ("c_", "s_"):
                    df.loc[tax_id, pre + "cumul"] = df.loc[tax_id, pre + "count"] + sum([
                        df.loc[child, pre + "cumul"]
                        for child in self.tax_id_hierarchy[tax_id]
                    ])
                    df.loc[tax_id, pre + "perc"] = 100 * df.loc[tax_id, pre + "cumul"] / total_counts

            # Format percentage signs.
            df.loc[:, "c_perc"] = round(df["c_perc"], 3).map(str) + "%"
            df.loc[:, "s_perc"] = round(df["s_perc"], 3).map(str) + "%"

            # Employ cutoff.
            cutoff = max((self.cutoff, (total_counts / 1e6) * self.cpm))
            df = df[(df['c_cumul'] > cutoff) | (df['s_cumul'] > cutoff) | (df.index == 'unclassified')]

            df.to_csv(os.path.join(tax_dir, sample_name + ".csv"), sep="\t", header=True)

            #
            # Map raw read output to taxonomy.
            #

            raw_counts_file = os.path.join(raw_counts_dir, sample_name + ".csv")
            raw_read_data = []

            with open(raw_counts_file, 'r') as f:
                for line in f:
                    read_data = line.strip().split('\t')[:-1]

                    if read_data[0] == 'U':
                        read_data[2] = 'NA'
                    else:
                        phy_id = int(read_data[2].lstrip('p'))
                        phylogeny_node = self.phylogeny_index.pool[phy_id].name
                        read_data[2] = phylogeny_to_taxonomy[phylogeny_node]

                    raw_read_data.append('\t'.join(read_data))

            with open(os.path.join(raw_output_dir, sample_name + '.csv'), 'w') as f:
                f.write('\n'.join(raw_read_data))

    def map_phylogeny_node(self, node_name, name_to_lineage, taxon_to_rank):
        # Map from (node_name) --> (taxid, rank, taxonomic lineage)
        #
        leaves = self.get_leaves(node_name)
        lineages = [name_to_lineage[leaf] for leaf in leaves]
        common_lineage = self.tuple_intersect(*lineages)

        if not common_lineage:
            phy_tax_id, phy_rank = "1", "root"
        else:
            phy_tax_id, phy_rank = taxon_to_rank[common_lineage[-1]]

        tax_id_map = {}
        pool_index = 1
        for i in range(len(common_lineage)):
            parent_tax_id, parent_rank = taxon_to_rank[common_lineage[i]]
            tax_id_map[parent_tax_id] = [parent_rank, " ".join(common_lineage[:i + 1])]

            if i == 0:
                # Connect to root.
                self.tax_id_hierarchy["1"].add(parent_tax_id)

            # Set tax_id in pool.
            try:
                pool_index = self.tax_id_pool.index(parent_tax_id) + 1
            except ValueError:
                self.tax_id_pool.insert(pool_index, parent_tax_id)
                pool_index += 1

            # Set taxa in hierarchy.
            if i + 1 < len(common_lineage):
                child_tax_id, child_rank = taxon_to_rank[common_lineage[i + 1]]
                if parent_tax_id in self.tax_id_hierarchy:
                    self.tax_id_hierarchy[parent_tax_id].add(child_tax_id)
                else:
                    self.tax_id_hierarchy[parent_tax_id] = {child_tax_id}
            elif parent_tax_id not in self.tax_id_hierarchy:
                self.tax_id_hierarchy[parent_tax_id] = set()

        return phy_tax_id, tax_id_map

    def get_leaves(self, node_name):
        return self.phylogeny_index.get_child_leaves(node_name)

    @staticmethod
    def tuple_intersect(*tuples):
        intersection = ()

        length = min([len(tup) for tup in tuples])
        for i in range(length):
            item = tuples[0][i]

            for j in range(1, len(tuples)):
                if tuples[j][i] != item:
                    return intersection

            else:
                intersection += (item,)

        return intersection

    def draw_results(self):
        # Draw classified tree.
        self.phylogeny_index.draw_results(
            os.path.join(self.in_dir, "classified_counts.csv"),
            os.path.join(self.out_dir, "phylotree_classified.pdf"),
            skiprows=[1],
            groups=self.groups,
            cutoff=self.cutoff,
            cpm=self.cpm,
            colour_list=self.colour_list,
            name_to_taxon=self.name_taxa,
            use_phyla=self.phyla,
            keep_zeros=self.keep_zeros,
            use_node_names=self.use_node_names
        )

        # Draw unclassified tree.
        self.phylogeny_index.draw_results(
            os.path.join(self.in_dir, "splits_counts.csv"),
            os.path.join(self.out_dir, "phylotree_splits.pdf"),
            groups=self.groups,
            cutoff=self.cutoff,
            cpm=self.cpm,
            colour_list=self.colour_list,
            name_to_taxon=self.name_taxa,
            use_phyla=self.phyla,
            keep_zeros=self.keep_zeros,
            use_node_names=self.use_node_names
        )
