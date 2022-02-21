#!/usr/bin/env python

import random

import numpy as np

from expam.c.extract import reverse_complement_combine, \
    extract_kmers_from_read, reverse_complement


#random.seed(100)
alpha = ['A', 'C', 'G', 'T']


def rev_comp_check():
    k = 16
    max_kmer = (2**(2*k))-1

    kmer = np.ndarray((2, ), dtype=np.uint64)
    rev_kmer = np.zeros((1, ), dtype=np.uint64)

    kmer[1] = np.random.randint(low=0, high=max_kmer, dtype=np.uint64)
    reverse_complement(kmer, rev_kmer, k)

    format_bits = '0%db' % k

    print(format(kmer[1], "064b"))
    print(format(rev_kmer[0], "064b"))


def main():
    n = 250
    k = 31

    random_string = "".join(random.choice(alpha) for _ in range(n)).encode('utf8')

    res = reverse_complement_combine(random_string, random_string)
    kmers = extract_kmers_from_read(res, k)

    rc_kmers = kmers.copy()
    reverse_complement(k, rc_kmers)

    rc_kmers = rc_kmers[::-1]

    return np.all(kmers == rc_kmers)


if __name__ == '__main__':
    rev_comp_check()
