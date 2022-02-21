import random

import expam.c.extract
import expam.c.map
import numpy as np

import extracting_kmers

READ_LENGTH = 150
ALPHA = {"A", "C", "G", "T"}
ALPHA_MAP = {
    0: b"A",
    1: b"C",
    2: b"G",
    3: b"T"
}

BAD_SEQUENCE = "data/GCF_000218445.fna"
BAD_READS = [
    b'TAAAGCTTTCTTACGTTCCTGAAGTCCGAGAAATGCCCACTCTGTCATATCTGTCCATTTATGGATCGTACGGTCATACACCTGAACGTATCCGATCCTTGTCGGATGCATCCGCACGCTGTTCGACTGGTATTCTTTTCCGGTGTACGG',
    b'ATCCATAACTTTCCCGTCTCTGTTATCTCTTTATTCAGACTGTCCGAATGCCTTCGGCATCAGAATTATTTTTCCTGATTTCTTCTTCTGAATCCGAATACGAACGCTGCCGCTGCAAGGAGTAATGTTTTACGGTCAACGAAAACATGC',
    b'CACCTGGATCTGATTGGAATGGGTGCGATGATCATTGCTGAGATTGGCGACTTCAACCGTTTTAGTTCCCCTGACAAGATACTTGCCTATGCAGGACTTTCTCCCTCTACTTACCAATCAGGTCAGCTTGACGGCTGTTATTCTCACATG',
    b'TCTTCCCAACGGATAGCTTTTGTTTCTCCAACACGGATAAAGAGATAGAATTAATGAGTCTCTCACGAATTTAGTGTGAAAGCCGAACCTCATCTTGTATCAATGTAGTTCTTAAAAATTTTTCTCTTTAAAGAAATTTTCAAGGGTTGT',
    b'ACGATGAGAACAATATTTGCAGAATACAATCCACAGTGCAACAGCATTGATGTATATACTAATACGGGCTATATACTCCGTATTGATTGTTGGGAAGCAGAAAAAAAATTTAAGAACCACACCTGGATCTGATTGGAATGGGTGCGATGA',
    b'ATCAGAATTATTTTTCCTGATTTCTTCTTCTGAATCCGAATACGAACGCTGCCGCTGCAAGGAGTAATGTTTTACGGTCAACGAAAACATGCCCCTTCGGAGTCTGCTGTATATCGGACGCGCTTATGAAAGGCTTGTACCGCCGAGAAG',
    b'CTGATTAACAATCAAGCAAAGAGACTGTGGCTGGACACTTTTGTAGGTTTGTCAAACATATGTTACACCCCGATCAATAAATTCCCTCAGAACAGTTTGTGGTGATTTCCACCCCAGAGGTCTCATGGGAAACTGGTTATAATCCCTACG',
    b'ATATATCATAACCATACCTGTAGAAAAGGTAAAGAATGGTTATGATGGAAATAGCATACAGACCTTGGAGAGGGTCTATAAATACTACTACTTTATAATACGAAAAGATGCGAGCGTAACAAAATCGGTATCGCGTGGAAAATTCATTTT',
    b'ACTACTGGCTCATTCTTTGTATTATTGATTACTGACTCGCTTCCTGTTCCAGCACTGTGTTCAGTGGACGTTATCTGTCAGGAACAACACGGAAGCCTACAGAAATCCCGTAGTTCCTGAGAATATTTTAAAAATATAATCAGAAGAGTA',
    b'TTGAAAAACGAAAACATCGCCCCTTCGGAGTCTGCTGTATATCGGACGCGCTTATGAAAGGCTTGTACCGCCGAGAAGCCGTTATAAGAAGAAAATCGTTTCCTTCCGACACCGGAATTTTACACATTCTATAACGGAAAAGAAAAGTGG',
    b'ATATACTAATACGGGCTATATACTCCGTATTGATTGTTGGGAAGCAGAAAAAAAATTTAAGAACCACACCTGGATCTGATTGGAATGGGTGCGATGATCATTGCTGAGATTGGCGACTTCAACCGTTTTAGTTCCCCTGACAAGATACTT',
    b'GCTTAGATGCCTTCAGCGTTTATCCCTTCCCGGCTTGGCTACTCTGCTATGCACTTGGTAATACAACAGATACACCAGCGGCCAGTCCACCCCGGTCCTCTCGTACTAAGGGCAGCTCCTCTCAAATATCCTGAATGATTTACATTTACC',
    b'ACAGATAAGGAAGAAATCTTTTGAAAGAGGTCATGGAAACCATTCCGCATTGTTTGGAATAGTAAATGCGTTTAAAAAAGCACTTCTCCTTCTGAGTTTGTGATAAAACAATCGTGTTTGTCTTTTGCAACATCAATACCAATATAAATC',
    b'GGTACATCTCCTTTCGGGTGGTTGGTTTATAGTTTCTAGGCAACTCTATTTTACCACAATCCGTTGAGGAGATGTTTTGTTTTATAAAAGTTAAAATGCCGTATTTATGCGGCTTGTGGCGTTTCGCAAACACCTCTTTGCGGGAACACA',
    b'ACGACTGTAATTGTTCTTCGTTGTGATAAGATGTGTTTGAGTAGGACTGTAAGTTTACTCCAGACATCATCATAAGAGCAATCGTATGGGCATCCACTTTATCGGTTTTTGTCTGTCTAAGGCTAAGACTTTTACGAAACAGATTGGTAT',
    b'CCCCATGGAAAATATTCCTGCCGACTACTGGCTCATTCTTTGTATTATTGATTACTGACTCGCTTCCTGTTCCAGCACTGTGTTCAGTGGACGTTATCTGTCAGGAACAACACGGAAGCCTACAGAAATCCCGTAGTTCCTGAGAATATT',
    b'AGTCAATATCTTCCCAACGGATAGCTTTTGTTTCTCCAACACGGATAAAGAGATAGAATTAATGAGTCTCTCACGAATTTAGTGTGAAAGCCGAACCTCATCTTGTATCAATGTAGTTCTTAAAAATTTTTCTCTTTAAAGAAATTTTCA',
    b'TAAGAAAGCAGGGAAGAGAAGAAGGGCAGAAGAAGGGAAGAGAAGAAGGACGAATAGAAGAAAAAAGTGCTCTCATCCGGAAAAAACTTGAAAAAGGAAAAACAATTTCCGAAATTGCAGATGATCTTGAAGATACAGAAGAAAATATTG',
    b'TGGTTTATAGTTTCTAGGCAACTCTATTTTACCACAATCCGTTGAGGAGATGTTTTGTTTTATAAAAGTTAAAATGCCGTATTTATGCGGCTTGTGGCGTTTCGCAAACACCTCTTTGCGGGAACACAGGACGAGGTGAATCAGATGGCA',
    b'ATTTTGGGTTTTACGACCGGCAATACAGAAACTGTCTTTAAAATGGTTGATTGATTTTTCAACGTTTTACACGGATTTTATAGGTCTCATCCCATTCTACGGAACCACGTTTCCACACCGGGATAGGCACGCAGATTCTTTCACCTGTTT'
]


def get_sequence():
    # Get main sequence.
    sequence = ""
    with open(BAD_SEQUENCE, "r") as f:
        for line in f:
            if ">" not in line:
                sequence += line.strip()

            else:
                sequence += "N"

    return bytes(sequence, "utf8")


def non_acgt(sequence):
    return [
        1 if nt in ALPHA else 0
        for nt in sequence.decode("utf8")
    ]


def check_increasing(arr):
    return np.all(arr[1:, :] > arr[:-1, :], axis=0).tolist()[0]


def make_read():
    num_n = random.randint(0, 10)
    sequence = [
        ALPHA_MAP[random.randint(0, 3)]
        for _ in range(READ_LENGTH - num_n)
    ]

    for _ in range(num_n):
        insertion_index = random.randint(0, len(sequence))
        sequence.insert(insertion_index, b"N")

    return b"".join(sequence)


def count_kmers(sequence, k):
    read_mask = non_acgt(sequence)

    for i in range(len(read_mask)):
        if read_mask[i] == 0:
            for j in range(k):
                if i - j >= 0:
                    read_mask[i - j] = 0

    for i in range(1, k):
        read_mask[-i] = 0

    return sum(read_mask)


def check_read_exists(sequence):
    for i, read in enumerate(BAD_READS):
        if read in sequence:
            print("(%d / %d)\t\tSucceeded!"
                  % (i + 1, len(BAD_READS)))

        else:
            print("(%d / %d)\t\tFailed!"
                  % (i + 1, len(BAD_READS)))


def check_reads(k, reference_kmers, values, null_value):
    bad_indices = []

    for i, read in enumerate(BAD_READS):
        # Extract kmers.
        kmers = expam.c.extract.extract_kmers_from_read(read, k)

        # Search against references.
        mapped_kmers = expam.c.map.binarysearch_get(
            kmers,
            reference_kmers,
            values,
            null_value
        )

        if null_value in mapped_kmers:
            print("(%d / %d)\t\tFailed!\t\t(Rank %d.)"
                  % (i + 1, len(BAD_READS), int(np.sum(mapped_kmers == null_value))))

            bad_indices.append(i)

        else:
            print("(%d / %d)\t\tSucceeded!"
                  % (i + 1, len(BAD_READS)))

    return bad_indices


def check_extractions(k, indices):
    for i, index in enumerate(indices):
        read = BAD_READS[index]

        base, _ = extracting_kmers.get_kmers(read, k)
        this_reference = base[:, 1:]
        this_kmers = expam.c.extract.extract_kmers_from_read(read, k)
        values = np.ones(this_kmers.shape[0], np.uint16)
        null_value = 0

        mapped_kmers = expam.c.map.binarysearch_get(
            this_kmers,
            this_reference,
            values,
            null_value
        )

        if null_value in mapped_kmers:
            print("(%d / %d)\t\tFailed!\t\t(Rank %d.)"
                  % (i + 1, len(indices), int(np.sum(mapped_kmers == null_value))))

        else:
            print("(%d / %d)\t\tSucceeded!"
                  % (i + 1, len(indices)))


def check_genome_extraction(k):
    sequence = get_sequence()

    # Get raw list of kmers.
    raw_kmers = expam.c.extract.extract_kmers_from_read(sequence, k)

    # Get ordered list of kmers using numpy.
    base, _ = extracting_kmers.get_kmers(sequence, k)
    expam_kmers = base[:, 1:]

    # Get ordered list of kmers from numpy.
    numpy_kmers = np.unique(raw_kmers)
    numpy_kmers = numpy_kmers.reshape((numpy_kmers.shape[0], 1))

    print("Numpy agrees with exPAM? %s" % str(np.array_equal(expam_kmers, numpy_kmers)))


def check_map_function(k, reference_kmers):
    all_kmers = reference_kmers.ravel()
    print("Reference is increasing? %s\n" % str(check_increasing(reference_kmers)))

    for i in range(len(BAD_READS)):
        read = BAD_READS[i]
        read_kmers = expam.c.extract.extract_kmers_from_read(read, k)

        search_indices = np.searchsorted(
            all_kmers,
            read_kmers
        )

        comparisons = np.take(all_kmers, search_indices, mode="clip") == read_kmers
        is_equal = np.all(comparisons)

        if not is_equal:
            print("(%d / %d)\t\tFailed\t\t(difference of %d / %d)!"
                  % (i + 1, len(BAD_READS), np.sum(comparisons), read_kmers.shape[0]))

        else:
            print("(%d / %d)\t\tSucceeded!"
                  % (i + 1, len(BAD_READS)))


def main():
    k = 31

    print("Getting sequence...")
    sequence = get_sequence()

    print("Getting reference kmers...")
    base, _ = extracting_kmers.get_kmers(sequence, k)
    reference_kmers = base[:, 1:]
    values = np.ones(reference_kmers.shape[0], dtype=np.uint16)
    null_value = 0

    print("Checking reference kmers...")
    print("\tStrictly increasing? %s" % str(check_increasing(reference_kmers)))

    print("\n*** Check 0: Reads exist in sequence ***")
    check_read_exists(sequence)

    print("\n*** Check 1: Compare to sorted references ***")
    # Get list of failed reads.
    failed_read_indices = check_reads(k, reference_kmers, values, null_value)

    print("\n*** Check 2: Compare extraction methods ***")
    # Compare extraction methods for these indices.
    check_extractions(k, failed_read_indices)

    print("\n*** Check 3: Check extraction of full genome ***")
    # Check my extraction agrees with numpy.
    check_genome_extraction(k)

    print("\n*** Check 4: Check binarysearch_get function ***")
    # Check my searching against numpy searching.
    check_map_function(k, reference_kmers)


if __name__ == '__main__':
    main()
