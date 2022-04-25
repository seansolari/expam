#cython: infer_types=True
#cython: boundscheck=False
#cython: wraparound=False
#cython: nonecheck=False
#cython: cdivision=True

cimport cython
from libc.math cimport ceil, floor
import numpy as np
cimport numpy as np

ctypedef np.uint8_t uint8_t
ctypedef np.uint16_t uint16_t
ctypedef np.uint64_t uint64_t
ctypedef np.int64_t int64_t


"""
Binary search implementation.

"""

cpdef binarysearch(uint64_t[:] target, uint64_t[:, :] kmers):
    """Binary search for target within array of kmers. Will return closest
    element in list if target is not found.

    :param target: Element to search for.
    :param kmers: Array to be searched within.
    :return: Whether element was found, and the corresponding index.
    :rtype: (int, int)

    """
    cdef:
        long long int i=0
        long long int L=0, U=kmers.shape[0]-1
        int j=0, W=kmers.shape[1]

    while True:
        if L > U:
            # Return the index with a fail signal.
            return 0, i
        i = <long long int>floor((U + L) / 2)
        for j in range(W):
            if kmers[i, j] > target[j]:
                U = i - 1
                break
            elif kmers[i, j] < target[j]:
                L = i + 1
                break
        else:
            return 1, i


def binarysearch_put(uint64_t[:, :] kmers, uint64_t[:, :] keys, uint16_t[:] values, uint16_t[:, :] lca_matrix, int new_id):
    """For all the elements in kmers, each coming from a node in the phylogeny
    with ID new_id, update their LCA values to incorporate this new_id. All
    possible LCA values have been computed in the LCA matrix.

    This means that every value in kmers should be in keys, as a binary search
    is executed to find the value of each kmer in keys.

    :param kmers: New items to be mapped.
    :param keys: (Sorted) Database keys to be searched against.
    :param values: Current LCA values for keys.
    :param lca_matrix: LCAs for each pair of nodes in the phylogeny.
    :param int new_id: ID of node containing kmers.
    :return: None (in-place)

    """
    cdef:
        long long int i, ind, M=kmers.shape[0]
        bint isin
        uint16_t current_id

    for i in range(M):
        # Find value in list.
        isin, ind = binarysearch(
            kmers[i, :],
            keys
        )

        # Find LCA of new and current value.
        current_id = values[ind]
        if current_id == new_id:
            pass
        elif current_id == 0:
            values[ind] = new_id
        elif current_id < new_id:
            values[ind] = lca_matrix[new_id, current_id]
        elif current_id > new_id:
            values[ind] = lca_matrix[current_id, new_id]


cpdef binarysearch_get(uint64_t[:, :] kmers, uint64_t[:, :] keys, uint16_t[:] values, uint16_t nullvalue):
    """Get the value of each element in kmers from keys. If the element is not
    found, nullvalue is used instead.

    :param kmers: Array of kmers to be searched for in keys.
    :param keys: Sorted database keys.
    :param values: Corresponding LCA values for keys.
    :param nullvalue: Value to be used if searched item is not found.
    :return: Array of values.
    :rtype: np.ndarray

    """
    cdef:
        long long int i, j, ind
        bint isin
        long long int M=kmers.shape[0]

    results = np.empty(M, dtype=np.uint16)
    cdef uint16_t[:] results_view = results
    for i in range(M):
        isin, ind = binarysearch(
            kmers[i, :],
            keys
        )
        if isin == 0:
            results_view[i] = nullvalue
        else:
            results_view[i] = values[ind]
    return results


"""
Insertion sort implementation.

"""

cdef unique(uint16_t[:] nodes, uint16_t[:] counts):
    insertion_sort(nodes, counts)
    nodes, counts = remove_duplicates(nodes, counts)
    return nodes, counts

cdef insertion_sort(uint16_t[:] nodes, uint16_t[:] counts):
    cdef int i = 1, j

    while i < nodes.shape[0]:
        j = i
        while j > 0 and nodes[j - 1] > nodes[j]:
            swap(nodes, j - 1, j)
            swap(counts, j - 1, j)
            j -= 1

        i += 1

cdef swap(uint16_t[:] arr, int i, int j):
    cdef uint16_t temp = arr[i]

    arr[i] = arr[j]
    arr[j] = temp

cdef remove_duplicates(uint16_t[:] nodes, uint16_t[:] counts):
    cdef int i = 0, j = 1

    while j < nodes.shape[0]:
        if nodes[i] == nodes[j]:
            counts[i] += counts[j]

        else:
            i += 1
            nodes[i] = nodes[j]
            counts[i] = counts[j]

        j += 1

    return nodes[:i + 1], counts[:i + 1]


"""
Classification routines.

"""

cpdef classify(uint64_t[:, :] kmers, uint64_t[:, :] keys, uint16_t[:] values, uint16_t[:, :] lca_matrix,
               uint16_t null_value, float alpha):
    cdef:
        uint16_t[:] nodes, counts
        uint16_t[:] raw_mapped_nodes
        uint16_t cls
        int candidates, cls_code = 0         # Innocent until proven guilty for class code.
        int not_null, required_counts, terminals
        str comp_string = ""

    # Map kmers using database.
    raw_mapped_nodes = binarysearch_get(kmers[:, 1:], keys, values, null_value)
    not_null = remove_null(raw_mapped_nodes, null_value)

    # Catch case of all unclassified.
    if not_null == 0:
        return -1, null_value, comp_string

    elif not_null < raw_mapped_nodes.shape[0]:
        cls_code = 1

    raw_mapped_nodes = raw_mapped_nodes[:not_null]

    # Convert mapped read to compressed string format.
    nodes, counts = kmers_to_compressed_form(raw_mapped_nodes)
    comp_string = compressed_to_string(nodes, counts)

    # Classification from lineage(s), with cutoff from alpha.
    required_counts = <int>ceil(alpha * kmers.shape[0])
    nodes, counts = unique(nodes, counts)

    # Get list of terminals.
    terminals = accumulate_candidates(nodes, counts, lca_matrix)

    # Find LCA of valid candidates.
    candidates, cls = choose_lca_candidate(nodes, counts, lca_matrix, terminals, required_counts)

    if candidates == 0:
        cls_code = 1
        cls = 1

    elif candidates > 1:
        cls_code = 1

    return cls_code, cls, comp_string


cdef int remove_null(uint16_t[:] raw_mapped_nodes, uint16_t null_value):
    cdef int j = 0, n = raw_mapped_nodes.shape[0], i

    for i in range(n):
        if raw_mapped_nodes[i] != null_value:
            raw_mapped_nodes[j] = raw_mapped_nodes[i]
            j += 1

    return j


cdef kmers_to_compressed_form(uint16_t[:] raw_mapped_nodes):
    """
    Return two arrays:
        1. Array of phylogeny nodes that consecutive kmers mapped to.
        2. The number of consecutive kmers that mapped to this node.

    """
    cdef:
        uint16_t[:] nodes, counts
        int n_kmers = raw_mapped_nodes.shape[0]
        int i, j = 0

    nodes = np.ndarray(n_kmers, dtype=np.uint16)
    counts = np.ones(n_kmers, dtype=np.uint16)
    
    nodes[0] = raw_mapped_nodes[0]

    for i in range(1, n_kmers):
        if raw_mapped_nodes[i] == nodes[j]:
            counts[j] += 1
        else:
            j += 1
            nodes[j] = raw_mapped_nodes[i]

    return nodes[:j + 1], counts[:j + 1]


def compressed_to_string(nodes, counts):
    return " ".join("p%d:%d" % (nodes[i], counts[i]) for i in range(nodes.shape[0]))


cdef int accumulate_candidates(uint16_t[:] nodes, uint16_t[:] counts, uint16_t[:, :] lca_matrix):
    cdef:
        int i, terminals = 0, n_nodes = nodes.shape[0]
        bint new_terminal

    for i in range(n_nodes - 1, -1, -1):
        new_terminal = check_terminals(nodes, counts, lca_matrix, i, terminals)

        if new_terminal:
            terminals = add_terminal(nodes, counts, i, terminals)

    return terminals


cdef bint check_terminals(uint16_t[:] nodes, uint16_t[:] counts, uint16_t[:, :] lca_matrix, int i, int terminals):
    cdef:
        int j, n = nodes.shape[0] - 1

    for j in range(n, n - terminals, -1):
        if lca_matrix[nodes[j], nodes[i]] == nodes[i]:          # Remember that j > i (Indexing LT lca matrix).
            counts[j] += counts[i]
            return 0

    return 1


cdef int add_terminal(uint16_t[:] nodes, uint16_t[:] counts, int i, int terminals):
    cdef:
        int minus_t = (nodes.shape[0] - 1) - terminals

    nodes[minus_t] = nodes[i]
    counts[minus_t] = counts[i]

    return terminals + 1


cdef choose_lca_candidate(uint16_t[:] nodes, uint16_t[:] counts, uint16_t[:, :] lca_matrix, int terminals,
                                   int required_counts):
    cdef:
        int n = nodes.shape[0] - 1, start
        int candidates = 0
        uint16_t lca = 1

    # Base case.
    start = n
    while start > n - terminals:
        if counts[start] >= required_counts:
            candidates = 1
            lca = nodes[start]

            break

        start -= 1

    # Take LCA of remaining valid candidates.
    for i in range(start - 1, n - terminals, -1):
        if counts[i] >= required_counts:
            candidates += 1

            if lca < nodes[i]:
                lca = lca_matrix[nodes[i], lca]
            else:   # Note that lca is never equal to any terminal in `nodes`.
                lca = lca_matrix[lca, nodes[i]]

    return candidates, lca
