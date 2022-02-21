#!/home/ssol0002/expam/venv/bin/python

import numpy as np
import expam.c.extract
from expam.run import Timer

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

def do_kmers(f_url, k, kind):
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

    sequence = get_sequence_data(furl)
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

    view.sort(order=['k%d' % chunk for chunk in range(chunks)], kind=kind)
    new_size = expam.c.extract.remove_duplicates(base)
    resize_base(new_size)

def mean(arr):
    return sum(arr) / len(arr)

def std(arr):
    av = mean(arr)
    return sum([abs(v - av) for v in arr]) / len(arr)

def main(f_url, k):
    data = {
        "quicksort": [],
        "mergesort": [],
        "heapsort": [],
        "stable": [],
    }

    for sort_kind, arr in data.items():
        print(sort_kind, "beginning...")

        for _ in range(5):
            with Timer() as t:
                do_kmers(f_url, k, sort_kind)

            arr.append(float(t))

    for sort_kind, arr in data.items():
        print("%s\t%fpm%f s" % (sort_kind, mean(arr), std(arr)))

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
