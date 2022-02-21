#!/usr/bin/env python3
import math
import random

from expam.c.extract import get_raw_kmers
import expam.c.map
from expam.c.map import classify, new_classify
from expam.run import Timer, timeit
import expam.tree
import numpy as np


ALPHA = ['A', 'C', 'G', 'T']
NULL_VALUE = 0

"""

"""


def format_mapped_kmers(arr):
    result = ""
    current_kmer, current_count = arr[0], 1

    i = 1
    while i < len(arr):
        if arr[i] == current_kmer:
            current_count += 1

        else:
            result += "p%d:%d " % (current_kmer, current_count)

            current_kmer = arr[i]
            current_count = 1

        i += 1

    result += "p%d:%d" % (current_kmer, current_count)

    return result


def old_process_read(read, k, read_array, keys, values, lca_matrix):
    with Timer() as t:
        read_length = len(read)
        n_kmers = get_raw_kmers(read, k, read_array)

        # Control for reads that don't have long enough good regions for kmers to be extracted.
        if n_kmers == 0:
            return "U", 0, read_length, ""

        # Map the kmers.
        cls, common_clade, mapped_kmers = classify(k, read_array[:n_kmers], keys,
                                                   values, lca_matrix, NULL_VALUE)
        kmer_string = format_mapped_kmers(mapped_kmers)

        if cls == -1:
            result = "U"

        elif cls == 1:
            result = "S"

        else:
            result = "C"

    return float(t)


def new_process_read(read, k, read_array, keys, values, lca_matrix):
    with Timer() as t:
        read_length = len(read)
        n_kmers = get_raw_kmers(read, k, read_array)

        # Control for reads that don't have long enough good regions for kmers to be extracted.
        if n_kmers == 0:
            return "U", 0, read_length, ""

        # Map the kmers.
        cls, common_clade, kmer_string = new_classify(read_array[:n_kmers], keys, values,
                                                      lca_matrix, NULL_VALUE, 0.5)

        if cls == -1:
            result = "U"

        elif cls == 1:
            result = "S"

        else:
            result = "C"

    return float(t)


def check_different_classification_methods():
    rng = np.random.default_rng(seed=100)

    n_trials = 1000

    k = 31
    n_leaves = 6
    n_nodes = n_leaves * (n_leaves - 1) / 2
    read_len = 150

    # Generate read of length 150.
    read = bytes(''.join(rng.choice(ALPHA) for _ in range(read_len)), encoding='utf8')

    # Create array for kmers.
    store = np.ndarray((read_len, 2), dtype=np.uint64)

    # Create key-value store.
    pre_keys = np.ndarray((read_len, 2), dtype=np.uint64)
    n_kmers = get_raw_kmers(read, k, pre_keys)
    kmers = pre_keys[:n_kmers, 1:]

    keys = kmers.astype(np.uint64)
    keys.sort(axis=0)
    values = rng.choice(np.arange(n_nodes, dtype=np.uint16), size=keys.shape[0])

    # Create LCA matrix.
    lca_matrix = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                           [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                           [0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                           [0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0],
                           [0, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0],
                           [0, 1, 2, 3, 4, 4, 0, 0, 0, 0, 0, 0],
                           [0, 1, 2, 3, 4, 4, 6, 0, 0, 0, 0, 0],
                           [0, 1, 2, 3, 4, 4, 6, 6, 0, 0, 0, 0],
                           [0, 1, 2, 3, 3, 3, 3, 3, 3, 0, 0, 0],
                           [0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0],
                           [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]], dtype=np.uint16)
    old_times = [
        old_process_read(read, k, store, keys, values, lca_matrix)
        for _ in range(n_trials)
    ]
    new_times = [
        new_process_read(read, k, store, keys, values, lca_matrix)
        for _ in range(n_trials)
    ]

    print("Old time taken: (%.2f +/- %.2f) x 1e-5" % (mean(old_times) * 1e5, sd(old_times) * 1e5))
    print("New time taken: (%.2f +/- %.2f) x 1e-5" % (mean(new_times)*1e5, sd(new_times)*1e5))


def mean(arr):
    return sum(arr) / len(arr)


def sd(arr):
    m = mean(arr)

    return math.sqrt(
        sum((x - m)**2 for x in arr) / len(arr)
    )


@timeit
def expam_unique(values, counts):
    expam_values, expam_counts = unique(values, counts)
    expam_values = np.array(expam_values, dtype=np.uint16)
    expam_counts = np.array(expam_counts, dtype=np.uint16)

    print(expam_values)
    print(expam_counts)


@timeit
def numpy_unique(values):
    numpy_values, numpy_counts = np.unique(values, return_counts=True)
    print(numpy_values)
    print(numpy_counts)


def check_insertion():
    rng = np.random.default_rng(seed=1)
    size = 10

    values_one = rng.choice(np.arange(10, dtype=np.uint16), size=size)
    values_two = np.copy(values_one)
    counts = np.ones(size, dtype=np.uint16)

    expam_unique(values_one, counts)
    print("\n...\n")
    numpy_unique(values_two)


"""

"""


def generate_random_tuple(a, b):
    return tuple(random.randint(0, 1) for _ in range(random.randint(a, b)))


def check_lca():
    child_one = generate_random_tuple(10, 100)
    child_two = generate_random_tuple(10, 100)

    lca = expam.tree.propose_lca(child_one, child_two)
    lca_length = len(lca)
    print("LCA length %d." % lca_length)

    if (lca_length > 0) and not (child_one[-lca_length:] == child_two[-lca_length:]):
        print("Failed!")

    if child_one[-lca_length-1] == child_two[-lca_length-1]:
        print("Failed!")

    print("Succeeded!")


def check_binary_search():
    n = random.randint(100, 1000)
    correct = 0

    keys = np.arange(n, dtype=np.uint64).reshape((n, 1))

    for i in range(n):
        isin, ind = expam.c.map.binarysearch(keys[i], keys)

        if (not isin) or (ind != i):
            print("Index %d failed!" % i)

        else:
            correct += 1

    else:
        print("binarysearch succeeded %d / %d!" % (correct, n))


def check_binary_get():
    n = random.randint(100, 1000)

    keys = np.arange(n, dtype=np.uint64).reshape((n, 1)) + 1
    values = np.arange(n, dtype=np.uint16) + 1
    nullvalue = 0

    results = expam.c.map.binarysearch_get(
        kmers=keys,
        keys=keys,
        values=values,
        nullvalue=nullvalue
    )

    if not np.all(results == values):
        print("binarysearch_get failed %d / %d!" % (np.sum(results != values), n))

    else:
        print("binarysearch_get succeeded %d / %d!" % (n, n))


"""

"""


def main():
    # check_binary_search()
    # check_binary_get()

    # check_lca()

    # check_insertion()
    check_different_classification_methods()


if __name__ == '__main__':
    main()
