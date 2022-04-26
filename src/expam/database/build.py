from ctypes import c_uint64
import math
from multiprocessing import Pipe, Value, shared_memory
import os
import subprocess
from typing import List, Mapping
import numpy as np
from expam.database import CHUNK_SIZE, TIMEOUT, UNION_RATIO, FileLocationConfig, expam_dtypes
from expam.database.config import load_database_config
from expam.process.genome import ExtractWorker
from expam.process.manager import ControlCenter, ExpamProcesses
from expam.process.piler import UnionWorker
from expam.sequences import COMP_PARSE, check_suffix
from expam.tree.tree import Index, propose_lca
from expam.utils import ls


def main(
    db_path: str, genome_paths: List[str], phylogeny_path: str, k: int,
    n=None, n_extract=None, n_union=None, pile_size=None,
):
    # Configure number of processes.
    if n_union is not None and n_extract is not None:
        n_union, n_extract = int(n_union), int(n_extract)
    elif n is not None:
        n_union, n_extract = distribute_processes(n)
    else:
        raise ValueError("Number of processes is can't be defined!")

    # Processing parameter details.
    num_cols = calculate_n_chunks(k)
    ranges = make_ranges(n_union, k)
    genome_paths, maxsize = sort_by_size(genome_paths)

    # Create shm_allocations for in-place kmer processing.
    shm_allocations = prepare_kmer_allocations(maxsize, num_cols, n_extract)

    # Import phylogeny.
    _, index = import_phylogeny(phylogeny_path)
    # Create LCA matrix.
    node_to_index = {node.name: i for i, node in enumerate(index.pool) if i > 0}
    lca_matrix = make_lca_matrix(index, node_to_index)

    # Prepare pipes between main and union workers.
    main_cons, union_cons = make_pipes(n_union)
    # Lock system to prevent concurrent access to HDF5 file.
    lock_value = Value(c_uint64, lock=True)

    # Multiprocessing configuration.
    config: FileLocationConfig = load_database_config(db_path)
    mp_config = {
        "name": "expam",    # Job system title.
        "phases": [         # Two-pass algorithm.
            "import",       # Import sequences and take union.
            "map"           # Map kmer to LCA.
        ],
        "layers": {         # Top layer of child processes - extract sequence & kmers.
            "class": ExtractWorker,
            "class_args": { # Positional arguments to init class.
                "k": k,
                "shm_kmer_params": shm_allocations,
                "kmer_ranges": ranges,
                "node_map": node_to_index,
                "out_dir": db_path,
                "logging_dir": config.logs,
            },
            "parent": ControlCenter,                        # Receive work from main process.
            "child": {
                "class": UnionWorker,       # Take union of sets of kmers.
                "class_args": { # Each either length n_union or 1.
                    "main_con": union_cons,
                    "lca_matrix": lca_matrix,
                    "pile_size": pile_size,
                    "logging_dir": config.logs,
                    "lock_value": lock_value,
                },
                "parent": ExtractWorker,    # Receive work from ExtractWorkers.
                "child": None,  # Bottom of process hierarchy.
                "n": n_union,   # Number of processes.
            },
            "n": n_extract,
        },
        "timeout": TIMEOUT,     # Time to block for when checking queues.
    }

    process_network = ExpamProcesses.from_method_dict(main_cons, db_path, config.logs, mp_config)

    # Load files for importing stage.
    for file_dir in genome_paths:
        process_network.load_job(
            phase="import",
            job=file_dir
        )

    process_network.run()

    # Save LCA matrix.
    np.save(config.lca_matrix, lca_matrix)

    # Concatenate accession ID files.
    concatenate_accession_ids(config.phylogeny)

def distribute_processes(n):
    n_union = max(math.floor(UNION_RATIO * n), 1)
    n_extract = max(n - n_union, 1)

    return n_union, n_extract

def make_ranges(n, k):
    # Calculate hash offset for the k value.
    maxvalue = maximum_hash(k)

    # Interpolate between offset and maxvalue.
    dx = maxvalue // n
    return {
        i: (i * dx, (i + 1) * dx) if i < n - 1
        else (i * dx, maxvalue)
        for i in range(n)
    }

def import_phylogeny(phylogeny_path):
    print("Importing phylogeny...")

    # Load the phylogeny into an Index.
    genome_names, index = Index.load_newick(phylogeny_path)
    return genome_names, index

def maximum_hash(k):
    if k < CHUNK_SIZE:
        return 4 ** k - 1
    else:
        return 4 ** CHUNK_SIZE - 1

def sort_by_size(dirs):
    _file_suffixes = COMP_PARSE
    _suffix_check = check_suffix

    def _gzip_size(file_dir):
        file_stats = subprocess.run(["gzip", "-l", file_dir], capture_output=True)
        file_meta = [entry for entry in file_stats.stdout.split(b' ') if entry]
        return file_meta[5]

    def _fna_size(file_dir):
        return os.stat(file_dir).st_size

    def _get_file_size(file_dir):
        for suffix in _file_suffixes:
            if _suffix_check(file_dir, suffix):
                return _gzip_size(file_dir)

        return _fna_size(file_dir)

    max_filename = max([len(f) for f in dirs])
    file_dt = np.dtype([
        ("filename", "S%d" % max_filename),
        ("filesize", np.int64)
    ])

    # Turn into structured array.
    file_array = np.array([
        (file_dir, _get_file_size(file_dir))
        for file_dir in dirs
    ], dtype=file_dt)

    # Sort by file size.
    file_array.sort(order="filesize")

    max_size = file_array["filesize"][-1]
    ordered_files = [
        filename.decode("utf8")
        for filename in file_array["filename"].tolist()
    ]

    # Assuming each char in the file is a base.
    return ordered_files, max_size

def prepare_kmer_allocations(rows, cols, n_processes):
    allocation_params = tuple()
    allocation_size = rows * cols * expam_dtypes.keys_dtype_size
    allocation_shape = (rows, cols)

    for _ in range(n_processes):
        # Create and 0 allocation.
        next_shm_allocation = shared_memory.SharedMemory(
            create=True,
            size=allocation_size
        )

        next_arr = np.ndarray(
            shape=allocation_shape,
            dtype=expam_dtypes.keys_dtype,
            buffer=next_shm_allocation.buf
        )
        next_arr[:] = 0

        # Pass on allocation details to children.
        allocation_params += (next_shm_allocation.name, allocation_shape)
        next_shm_allocation.close()

    return allocation_params

def calculate_n_chunks(k: int):
    return k // CHUNK_SIZE + (k % CHUNK_SIZE != 0)

def yield_coordinates(index: Index):
    n = len(index.pool)

    for i in range(1, n):
        for j in range(1, i):
            yield i, j, index.pool[i].coordinate, index.pool[j].coordinate

def compute_lca(i, j, coord_one, coord_two):
    return i, j, propose_lca(coord_one, coord_two)

def make_lca_matrix(index: Index, node_to_index: Mapping[str, str]):
    print("Creating LCA matrix...")

    def get_children(fixed_node, flexible_index_list):
        nonlocal index, matrix, node_to_index

        i = 0
        while i < len(flexible_index_list):
            for child in index.pool[flexible_index_list[i]].children:
                m = min(fixed_node, node_to_index[child])
                n = max(fixed_node, node_to_index[child])

                if matrix[n, m] == 0 and n != m:
                    flexible_index_list.append(node_to_index[child])

            i += 1

    n = len(index.pool)
    matrix = np.zeros((n, n), dtype=np.uint16)

    for i in range(1, n):
        subject = index.pool[i].coordinate
        subject_name = index.pool[i].name

        for j in range(1, i):
            if matrix[i, j] > 0:
                continue

            target = index.pool[j].coordinate
            target_name = index.pool[j].name

            lca_coordinate = propose_lca(subject, target)
            lca = index.coord(coordinate=lca_coordinate).name

            left_group = [i]
            right_group = [j]

            if lca not in [subject_name, target_name]:
                # Subject and target are on separate lineages.
                get_children(j, left_group)
                get_children(i, right_group)

            elif lca == subject_name:
                # Target is contained in the clade rooted at subject.
                get_children(i, right_group)

            elif lca == target_name:
                # Subject is contained in the clade rooted at target.
                get_children(j, left_group)

            for p in left_group:
                for q in right_group:
                    m = min(p, q)
                    n = max(p, q)

                    if n != m and matrix[n, m] == 0:
                        matrix[n, m] = node_to_index[lca]

    return matrix

def make_pipes(n):
    parent_cons, child_cons = [], []

    i = 0
    while i < n:
        parent_con, child_con = Pipe()

        parent_cons.append(parent_con)
        child_cons.append(child_con)

        i += 1

    return tuple(parent_cons), tuple(child_cons)

def concatenate_accession_ids(f_url):
    def get_data(file_url):
        with open(file_url, "r") as f:
            data = f.read()

        # Delete temporary data file.
        os.remove(file_url)

        return data

    def get_files(folder_url):
        return [
            file_url
            for file_url in ls(folder_url)
            if "accession_ids_" in file_url
        ]

    id_data = [
        get_data(os.path.join(f_url, file_name))
        for file_name in get_files(f_url)
    ]

    with open(os.path.join(f_url, "accession_ids.csv"), "w") as f:
        f.write("\n".join(id_data))



