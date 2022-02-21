import os

import expam.c.extract
import expam.c.map
import expam.run
import numpy as np

def check_range(kmers, lower, upper):
    if not ((kmers[0][0] >= lower) and (kmers[-1][0] <= upper)):
        print("\t\tFailed at range (%d, %d)!" % (lower, upper))

def send_kmers(kmers, ranges):
    target = np.ndarray(shape=(kmers.shape[1]), dtype=kmers.dtype)
    target[:] = ranges[0][0]

    total = 0
    n = len(ranges) - 1
    lower_ind = upper_ind = get_index(target, kmers, greater=False)

    for i, (lower_bound, upper_bound) in ranges.items():
        target[:] = upper_bound
        upper_ind = get_index(target, kmers)

        # Upper bound can only attain equality for maximum hash value.
        if i == n:
            while (upper_ind < kmers.shape[0]) and kmers_equal(target, kmers[upper_ind]):
                upper_ind += 1

        if lower_ind != upper_ind:
            check_range(kmers[lower_ind:upper_ind], lower_bound, upper_bound)

            total += upper_ind - lower_ind

        # Beginning next search.
        lower_ind = upper_ind

    return total

def get_index(target, kmers, greater=True):
    def _check_ind():
        nonlocal ind, kmers

        if (greater and ind == kmers.shape[0]) or ((not greater) and ind == 0):
            return True

        return False

    # Find closest element in list. Note that this element
    # may be less than the bound.
    _, ind = expam.c.map.binarysearch(target, kmers)

    if _check_ind():
        return int(ind)

    while (not _check_ind()) and kmers_compare(kmers[ind], target, greater):
        if greater:
            ind += 1

        else:
            ind -= 1

    return int(ind)

def kmers_equal(one, two):
    return np.all(one == two)

def kmers_compare(one, two, greater):
    for i in range(one.shape[0]):
        if one[i] < two[i]:
            return greater

        elif one[i] > two[i]:
            return not greater

        elif i+1 == one.shape[0]:
            return not greater

def main(seq_url, k, ranges):
    def resize_base(size):
        nonlocal base, chunks

        base = base[:size]
        update_view()

    def new_view():
        nonlocal base, d_type
        return base.ravel().view(dtype=d_type)

    def update_view():
        nonlocal base, d_type, view
        view = base.ravel().view(dtype=d_type)

    chunks = k // 32 + (k % 32 != 0)
    d_type = np.dtype(
        [
            ('index', np.uint64),
            *(
                ('k%d' % chunk, np.uint64)
                for chunk in range(chunks)
            )
        ]
    )

    print("Getting sequence data from %s." % seq_url)
    seq_data = ""

    with open(seq_url, "r") as f:
        for line in f:
            if ">" not in line:
                seq_data += line.strip()

            else:
                seq_data += "N"

    sequence = bytes(seq_data, encoding='utf8')

    print("\tGetting kmers...")
    base = np.ndarray(
        shape=(len(seq_data), chunks + 1),
        dtype=np.uint64,
        order="C",
    )
    view = new_view()

    new_size = expam.c.extract.get_raw_kmers(
        sequence=sequence,
        k=k,
        arr=base
    )
    resize_base(new_size)

    # Sort by kmers, which also sorts the indices
    view.sort(order=['k%d' % chunk for chunk in range(chunks)])
    new_size = expam.c.extract.remove_duplicates(base)
    resize_base(new_size)

    print("\tSending kmers...")
    n_sent = send_kmers(base[:, 1:], ranges)
    print("\t\t%d/%d sent!" % (n_sent, base.shape[0]))


if __name__ == '__main__':
    k = 31
    n = 15
    seqs = [
        "GCF_002240135.1_ASM224013v1_genomic.fna",
        "GCF_000021125.1_ASM2112v1_genomic.fna",
        "GCF_000055945.1_ASM5594v1_genomic.fna",
        "GCF_006874785.1_ASM687478v1_genomic.fna",
        "GCF_000513215.1_DB11_genomic.fna"
    ]
    seqs = [os.path.join("data", seq) for seq in seqs]

    ranges = expam.run.make_ranges(n, k)

    for seq in seqs:
        main(seq, k, ranges)
