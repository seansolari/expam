#!/usr/bin/env python3

"""
author: sean.solari@monash.edu
date: August, 2024

Given an Expam database and valid pointers (within the database configuration
file) to the reference sequences used to build the database, calculate:
    - The number of unique k-mers for each reference sequence, and
    - The average duplicity of the unique k-mers for each reference sequence.

IN:
    - Expam Database
OUT:
    - CSV table [genome, # k-mers, duplicity]

"""

from argparse import ArgumentParser
import gzip
import io
from multiprocessing import Pool, shared_memory
import os
from typing import Mapping
from expam.database import CHUNK_SIZE, expam_dtypes
from expam.database.config import JSONConfig, make_database_config
from expam.database.db import SharedDB
from expam.ext.kmers import get_raw_kmers
from expam.tree.tree import Index
import numpy as np
import pandas as pd

__version__ = (0, 1, 2)


COMPRESSION_EXTNS = ['.tar.gz', '.tar', '.gz']
DEFAULT_MODE = "rb"
DEFAULT_OPENER = open
COMP_PARSE = {
    ".tar.gz": {"mode": "rb", "opener": gzip.open},
    ".gz": {"mode": "rb", "opener": gzip.open}
}

def get_opener(file_name):
    for suffix in COMP_PARSE:
        if check_suffix(file_name, suffix):
            file_name = file_name.replace(suffix, '')
            mode, opener = COMP_PARSE[suffix]["mode"], COMP_PARSE[suffix]["opener"]
            break
    else:
        mode, opener = DEFAULT_MODE, DEFAULT_OPENER
    return file_name, (mode, opener)

def check_suffix(string, sfx):
    string = string.strip()
    if string[-len(sfx):] == sfx:
        return True
    return False

def format_name(name, remove_comp: bool = False):
    name = os.path.basename(name)
    # Must not be a digit. If it is, add 'ref' to the start.
    if name.isdigit():
        name = "ref" + name
    # Extract sequence name.
    if ".fna" in name:
        name = name.replace(".fna", "")
    elif ".faa" in name:
        name = name.replace(".faa", "")
    if remove_comp:
        for suffix in COMP_PARSE:
            if check_suffix(name, suffix):
                name = name.replace(suffix, '')
                break
    return name

def read_genome(genome_path: str, min_contig_size: int):
    # Check it is a valid file.
    if not os.path.isfile(genome_path):
        raise ValueError("Directory %s is not a file!" % genome_path)
    print("extracting sequences from %s..." % genome_path)
    # determine opening mode based on compression of sequence file
    file_name = format_name(genome_path)
    file_name, (mode, opener) = get_opener(file_name)
    # parse data
    read_stream = io.BytesIO()
    with opener(genome_path, mode=mode) as f:
        for line in f:
            if line.startswith(b'>'):
                read_stream.write(b'N')
            else:
                read_stream.write(line.strip())
    # finalise data
    read_str = read_stream.getvalue()
    print("finished reading sequences from %s" % genome_path)
    del read_stream
    return file_name, filter(lambda s: len(s) > min_contig_size, read_str.split(b'N'))

def calculate_n_chunks(k: int) -> int:
    return k // CHUNK_SIZE + (k % CHUNK_SIZE != 0)

def load_database(db_conf: JSONConfig):
    # get db params
    keys_shape, values_shape = tuple(db_conf["keys_shape"]), tuple(db_conf["values_shape"])
    if len(keys_shape) == 1:
        keys_shape = keys_shape + (1,)
    db_kmers = SharedDB(file_config.base, keys_shape, values_shape)
    return db_kmers

def attach_to_db(db: SharedDB):
    # load the database
    keys_shm = shared_memory.SharedMemory(
        name=db.keys_shm_name,
        create=False)
    keys = np.ndarray(
        shape=db.keys_meta["shape"],
        dtype=db.keys_meta["dtype"],
        buffer=keys_shm.buf)
    values_shm = shared_memory.SharedMemory(
        name=db.values_shm_name,
        create=False)
    values = np.ndarray(
        shape=db.values_meta["shape"],
        dtype=db.values_meta["dtype"],
        buffer=values_shm.buf)
    return keys_shm, keys, values_shm, values

def windowed_true(arr: np.ndarray, w: int):
    return np.lib.stride_tricks.sliding_window_view(arr, w).any(axis=1)

def calculate_kmer_duplicity(genome_path: str, genome_id: int, k: int, db_interface: SharedDB):
    assert calculate_n_chunks(k) == 1
    # attach to db interface
    keys_shm, keys, values_shm, values = attach_to_db(db_interface)
    # get k-mers that are unique to this genome
    unique_kmers = keys[values == genome_id, 0].ravel()
    unique_kmers_observed = 0
    genome_name, contigs = read_genome(genome_path, k)
    for sequence in contigs:
        # extract k-mers from sequence - note: sequence contains no 'N's
        sequence_length = len(sequence)
        arr = np.zeros(
            shape=(sequence_length, 2),   # account for index col
            dtype=expam_dtypes.keys_dtype,
            order="C")
        num_kmers = get_raw_kmers(sequence, k, arr)
        arr = arr[:num_kmers, 1].ravel()                        # ignore index column
        unique_kmers_observed += np.in1d(arr, unique_kmers)\
            .sum()
    # unattach from database
    keys_shm.close()
    values_shm.close()
    print("finished %s" % genome_name)
    return genome_name, unique_kmers.size, unique_kmers_observed / unique_kmers.size

def leaf_ids(tree: Index, name2file: Mapping[str, str]):
    return [
        (name2file[node.name], i)
        for i, node in enumerate(tree.pool)
        if i > 0 and node.type == "Leaf"
    ]

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--db", dest="database", type=str, required=True, help="Expam database")
    parser.add_argument("-p", "--parallel", dest="parallel", type=int, default=4, help="Number of processes to use")
    parser.add_argument("-o", "--out", dest="out", type=str, required=True, help="File where results are saved (csv format)")
    return parser.parse_args()

if __name__ == "__main__":
    print("Count Unique k-mers\n-------------------\n" + "\033[1m" + "\033[91m" + "Warning: " + "\033[0m" + "designed for small databases.\n")
    args = parse_args()
    file_config = make_database_config(args.database)
    db_conf = JSONConfig(file_config.conf)
    # get input genome paths
    k, _, phylogeny_path, genome_paths, _ = db_conf.get_build_params()
    # load expam database
    db_interface = load_database(db_conf)
    # get mapping from input genomes to values
    _, tree = Index.load_newick(phylogeny_path)
    name2file = {format_name(path, remove_comp=True): path for path in genome_paths}
    enumerated_genomes = leaf_ids(tree, name2file)
    # use paired or singular method
    print("searching for k-mer size k=%d" % k)
    fargs = (
        (genome_path, genome_id, k, db_interface)
        for genome_path, genome_id in enumerated_genomes
    )
    # calculate normalisation factors
    with Pool(processes=args.parallel) as p:
        data = p.starmap(calculate_kmer_duplicity, fargs)
    df = pd.DataFrame(data=data, columns=["genome_name", "unique_kmers", "kmer_dup"])
    df.sort_values(by="genome_name", inplace=True)
    df.to_csv(args.out, index=False)
    print("csv output written to %s" % args.out)
    # shut down SHM
    keys_shm, _, values_shm, _ = attach_to_db(db_interface)
    keys_shm.unlink()
    values_shm.unlink()
