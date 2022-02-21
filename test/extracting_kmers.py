import itertools

import numpy as np
import expam.c.extract

ALPHA = ["A", "C", "G", "T"]


def get_sequence_data(f_url):
    print("Getting sequence data.")
    seq_data = ""

    with open(f_url, "r") as f:
        for line in f:
            if ">" not in line:
                seq_data += line.strip()

            else:
                seq_data += "N"

    return bytes(seq_data, encoding='utf8')


def get_kmers(sequence, k):
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

    base = np.ndarray(
        shape=(len(sequence), chunks + 1),
        dtype=np.uint64,
        order="C",
    )
    view = new_view()

    # Create draft mask and get corresponding kmers.
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

    return base, view


def check_ordered(view):
    print("Final:")
    print("\t%d kmers." % view['index'].shape[0])
    print("\tIncreasing?", np.all(np.diff(view['k0']) > 0), sep=" ")
    _, s_counts = np.unique(view['k0'], return_counts=True)
    print("\tOrdered?", np.all(s_counts == 1), sep=" ")


def check_from_mask(sequence, base, view):
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

    print("Extracting mask.")
    mask = expam.c.extract.to_32bit(view['index'])

    print("Making new allocation.")
    big_next_base = np.ndarray(
        shape=(2 * mask.shape[0], chunks + 1),
        dtype=np.uint64,
        order="C",
    )
    next_base = big_next_base[:mask.shape[0]]
    next_view = next_base.ravel().view(dtype=d_type)

    print("Importing mask.")
    expam.c.extract.import_mask(
        data=mask,
        dest=next_view['index']
    )

    print("Extracting kmers from mask.")
    expam.c.extract.map_kmers(
        sequence=sequence,
        k=k,
        arr=next_base
    )

    print("Kmers equal?")
    print("\t", np.array_equal(base, next_base), sep="")
    print("\tIncreasing?", np.all(np.diff(next_view['k0']) > 0), sep=" ")
    _, s_counts = np.unique(next_view['k0'], return_counts=True)
    print("\tOrdered?", np.all(s_counts == 1), sep=" ")


def generate_ordered_sequence(l=5):
    return bytes("N".join([
        "".join(perm)
        for perm in itertools.product(ALPHA, repeat=l)
    ]), encoding="utf8")


def check_ordered_sequence():
    k = 6
    sequence = generate_ordered_sequence(k)

    base, view = get_kmers(sequence, k)

    validated = view["k0"] == np.arange(4 ** k, dtype=np.uint64)
    if not np.all(validated):
        print("kmer_extraction failed %d / %d!" % (np.sum(~validated), 4 ** k))

    else:
        print("kmer_extraction succeeded %d / %d!" % (4 ** k, 4 ** k))


def main(f_url, k):
    sequence = get_sequence_data(f_url)
    base, view = get_kmers(sequence, k)

    check_ordered(view)
    check_from_mask(sequence, base, view)
    check_ordered_sequence()


if __name__ == '__main__':
    from argparse import ArgumentParser
    from os.path import exists, isfile

    parser = ArgumentParser()
    parser.add_argument("-d", dest="dir")
    parser.add_argument("-k", dest="k")
    args = parser.parse_args()

    file_name, k = args.dir, int(args.k)

    if not (isfile(file_name) and exists(file_name)):
        print("File %s not found!" % file_name)

    else:
        main(file_name, k)
