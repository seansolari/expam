#cython: infer_types=True

cimport cython
from libc.stdlib cimport malloc, free
import numpy as np
cimport numpy as np

ctypedef np.uint64_t uint64_t
ctypedef np.int64_t int64_t

# Custom C scripts.
cdef extern from "mfil.c":
    cdef int64_t multidisjoint(int n_arrays, int n_columns, uint64_t **arr_ptrs, int64_t *sizes)
    cdef void k_way_merge(uint64_t **arr_ptrs, int64_t *indices, uint64_t n, uint64_t k, uint64_t chunks)

def filler(uint64_t k, uint64_t old_shm_size, list list_of_arrays):
    """Merge-sort all arrays contained in list_of_arrays into one array,
    namely the first element in this list of arrays. This means that the first
    list of list_of_arrays needs to be large enough to fit the distinct values
    of all other elements, however it can itself contain values.

    :param old_shm_size: The number of values contained in the first array.
    :param list_of_arrays: List containing arrays to be merged.
    :return: None, in-place on first array.
    :rtype: None

    """
    cdef:
        uint64_t n, n_arrays, chunks
        uint64_t **arr_ptrs
        int64_t[::1] indices
        uint64_t[:, ::1] temp_arr
        int i

    n = list_of_arrays[0].shape[0]
    n_arrays = len(list_of_arrays)
    chunks = list_of_arrays[0].shape[1]

    # Pointers on each of the arrays.
    arr_ptrs = <uint64_t **>malloc(k * sizeof(uint64_t *))

    if not arr_ptrs:
        raise MemoryError

    # Pointer to current index value.
    indices = np.zeros(k, dtype=np.int64)

    # Type cast each of the input lists and make pointer.
    for i in range(k):
        if i < n_arrays:
            temp_arr = list_of_arrays[i]  # Cast array.
            arr_ptrs[i] = &temp_arr[0, 0]  # Declare pointer on array.

            # Set the current index (bottom -> up merging).
            if i == 0:
                indices[i] = old_shm_size - 1
            else:
                indices[i] = temp_arr.shape[0] - 1

        # Add dummy indices until we can make a binary tree.
        else:
            arr_ptrs[i] = &temp_arr[0, 0]
            indices[i] = -1

    try:
        k_way_merge(arr_ptrs, &indices[0], n, k, chunks)

    finally:
        free(arr_ptrs)

def disjoint(list list_of_arrays):
    """Iteratively remove elements from the first array until the values left
    are distinct, with respect to the values contained in the other arrays.

    :param list_of_arrays: List of sorted arrays.
    :return: Number of distinct elements left in first array.
    :rtype: int

    """
    cdef:
        int n, chunks, i
        int64_t final_size
        int64_t[::1] sizes
        uint64_t **arr_ptrs
        uint64_t[:, ::1] temp_arr

    n = len(list_of_arrays)
    chunks = list_of_arrays[0].shape[1]
    final_size = 0  # Set default.

    sizes = np.zeros(n, dtype=np.int64)
    arr_ptrs = <uint64_t **>malloc(n * sizeof(uint64_t))

    if not arr_ptrs:
        raise MemoryError

    # Establish pointer on each array.
    for i in range(n):
        temp_arr = list_of_arrays[i]    # Cast array.
        arr_ptrs[i] = &temp_arr[0, 0]   # Declare pointer on array.

        # Establish array size.
        sizes[i] = list_of_arrays[i].shape[0]

    try:
        final_size = multidisjoint(n, chunks, arr_ptrs, &sizes[0])

    finally:
        # Free heap arrays.
        free(arr_ptrs)

    return final_size
