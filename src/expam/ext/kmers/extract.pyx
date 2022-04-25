#cython: infer_types=True
#cython: boundscheck=False
#cython: wraparound=False
#cython: nonecheck=False
#cython: cdivision=True

cimport cython
from libc.stdlib cimport malloc, free
cimport numpy as np
import numpy as np

ctypedef np.uint8_t uint8_t
ctypedef np.uint32_t uint32_t
ctypedef np.uint64_t uint64_t
ctypedef unsigned long long uint64
ctypedef long long int int64

cdef extern from "kmers.c":
    cdef uint64_t raw_kmers(const unsigned char *sequence, int k, uint64_t *arr)
    cdef void fill_kmers(const unsigned char *sequence, int k, uint64_t n, uint64_t *arr)

cdef extern from "jellyfish.c":
    cdef uint64_t word_reverse_complement(uint64_t word, int k)

cpdef void reverse_complement(uint64_t[::1] kmers, uint64_t[::1] r_kmers, int k):
    cdef:
        int j
        int m = kmers.shape[0]-1, b = 2 * ((k % 32) + (32 * (k % 32 == 0)))
        uint64_t r, mask = ((2**b) - 1)

    for j in range(m):
        if j == m-1:
            r_kmers[j] = word_reverse_complement(kmers[j+1], k)
        else:
            r_kmers[j] = word_reverse_complement(kmers[j+1], 32)

    if m > 1:
        # Swap order of chunks.
        for j in range(m // 2):
            r = r_kmers[m-(j+1)]
            r_kmers[m-(j+1)] = r_kmers[j]
            r_kmers[j] = r

        # Fix gaps in chunks.
        for j in range(m-1):
            r_kmers[j] <<= (64 - b)
            r = r_kmers[j+1] >> b
            r_kmers[j] |= r

        r_kmers[m-1] &= mask

cdef bint kmer_less(uint64_t[::1] a, uint64_t[::1] b):
    cdef int i, l = a.shape[0]

    for i in range(l):
        if a[i] < b[i+1]:
            return 1
        elif a[i] > b[i+1]:
            return 0

    return 0

cdef void canonical(uint64_t[:, ::1] kmers, int k):
    cdef:
        int i, l = kmers.shape[0], m = kmers.shape[1]-1
        uint64_t[::1] kmer, r_comp = np.ndarray(m, dtype=np.uint64)

    for i in range(l):
        kmer = kmers[i]
        reverse_complement(kmer, r_comp, k)

        if kmer_less(r_comp, kmer):
            kmers[i, 1:] = r_comp[:]

def get_raw_kmers(const unsigned char[::1] sequence, int k, uint64_t[:, ::1] arr):
    """Extract k-mers from sequence into the last (k//32+!!(k%32!=0)) columns of arr.
    The indices of these kmers (with respect to sequence) are stored in the
    first column of arr.

    :param sequence: Nucleotide sequence.
    :type sequence: bytes
    :param k: k value for analysis.
    :param arr: Array where kmers and their indices will be placed.
    :return: The number of rows of arr where kmers have been placed.
    :rtype: int

    """
    cdef uint64_t n_kmers

    n_kmers = raw_kmers(&sequence[0], k, &arr[0, 0])
    canonical(arr, k)

    return n_kmers

def bytes_intersect(const unsigned char[:] p, const unsigned char[:] q):
    cdef int l = min(len(p), len(q))

    for i in range(l):
        if p[i] != q[i]:
            return bytes(p[:i])
    else:
        return bytes(p[:l])

def reverse_complement_combine(const unsigned char[:] old, const unsigned char[:] new):
    cdef:
        int l, n, d

    l = len(old)
    n = len(new)
    d = n + l

    view = np.frombuffer(old, dtype='S1')
    dest = np.ndarray(d + 1, dtype='S1')

    for i in range(l):
        dest[i] = view[i]

    dest[l] = b'N'

    for i in range(n):
        if new[i] == b'A':
            dest[d-i] = b'T'
        elif new[i] == b'C':
            dest[d-i] = b'G'
        elif new[i] == b'G':
            dest[d-i] = b'C'
        elif new[i] == b'T':
            dest[d-i] = b'A'
        else:
            dest[d-i] = new[i]

    return dest.tobytes()

def init_index(const unsigned char[::1] sequence, int k, uint64_t[:, ::1] arr):
    """Equivalent to get_raw_kmers, however valid kmer indices are first found
    and then the kmers corresponding to these indices are then extracted
    independently. (Prefer to use get_raw_kmers.)

    :param sequence: Nucleotide sequence.
    :type sequence: bytes
    :param k: k value for analysis.
    :param arr: Array where kmers and their indices will be placed.
    :return: The number of rows of arr where kmers have been placed.
    :rtype: int

    """
    cdef uint64_t i

    i = start_mask(sequence, k, arr[:, 0])
    fill_kmers(&sequence[0], k, i, &arr[0, 0])

    canonical(arr, k)

    return i

def remove_duplicates(uint64_t[:, :] arr):
    """Remove (in-place) duplicate elements of a sorted arr. Equality
    is determined by values in the second column onwards.

    :param arr: Array of values.
    :return: Number of unique elements remaining.
    :rtype: int

    """
    cdef:
        uint64 arr_rows, arr_cols
        uint64 i, j, x

    arr_rows = arr.shape[0]
    arr_cols = arr.shape[1]

    x = 1

    for i in range(1, arr_rows):

        # Check last indices first as they're more likely to not agree.
        for j in range(arr_cols - 1, 0, -1):    # Don't look at index column.
            if arr[i, j] != arr[i - 1, j]:

                arr[x, :] = arr[i, :]
                x += 1

                break

    return x

def to_32bit(uint64_t[:] arr):
    return np.ascontiguousarray(arr, dtype=np.uint32)

def import_mask(uint32_t[::1] data, uint64_t[:] dest):
    move_32bit(data, dest)

def map_kmers(const unsigned char[::1] sequence, int k, uint64_t[:, ::1] arr):
    """Given an array whose first column is filled with indices of elements in
    some sequence, extract the kmers starting at each of these indices
    into the remaining columns of arr.

    :param sequence: Nucleotide sequence.
    :type sequence: bytes
    :param k: k value.
    :param arr: Array containing indices and space for kmers.
    :return: None, operations are done in-place.

    """
    cdef uint64_t n

    n = arr.shape[0]
    fill_kmers(&sequence[0], k, n, &arr[0, 0])
    canonical(arr, k)

def extract_kmers_from_read(const unsigned char[::1] read, int k):
    """Extract the kmers from a sequence into a new array.

    :param read: Short (~150bp) nucleotide sequence.
    :param int k: k value for analysis.
    :return: Array of kmers.
    :rtype: np.ndarray

    """
    cdef:
        uint64_t length, chunks, n_kmers
        uint64_t[:, ::1] kmers

    chunks = k // 32 + (k % 32 != 0)
    length = len(read)

    kmers = np.ndarray(
        shape=(length, chunks + 1),
        dtype=np.uint64
    )

    n_kmers = raw_kmers(&read[0], k, &kmers[0, 0])
    canonical(kmers, k)

    return kmers[:n_kmers]

cdef uint64_t start_mask(const unsigned char[::1] sequence, int k, uint64_t[:] index):
    cdef:
        uint64_t sequence_length
        uint64_t i, j, m
        char *valid_bits

    sequence_length = sequence.shape[0]
    i = j = 0

    valid_bits = <char *> malloc(sequence_length * sizeof(char))

    if not valid_bits:
        raise MemoryError

    try:
        # Each bit is a valid candidate index.
        for i in range(sequence_length):
            valid_bits[i] = 1

        for i in range(sequence_length):
            if sequence[i] not in b"ACGT":

                # Don't invalid bits using negative index.
                m = 0 if i + 1 - k < 0 else i + 1 - k

                for j in range(m, i+1):
                    valid_bits[j] = 0

        # Final edge.
        for j in range(k-1):
            valid_bits[-j] = 0

        # Convert bits into list of indices.
        j = 0
        for i in range(sequence_length):
            if valid_bits[i] == 1:
                index[j] = i
                j += 1

    finally:
        free(valid_bits)

    return j

cdef void move_32bit(uint32_t[::1] source, uint64_t[:] dest):
    cdef int i, length

    length = source.shape[0]

    for i in range(length):
        dest[i] = <uint64_t>source[i]
