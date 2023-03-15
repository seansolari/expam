#!/usr/bin/env python3
from argparse import ArgumentParser
import os

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--chunk", dest="chunk_size", type=int, required=True, help="Max batch size in Gb.")
    parser.add_argument("-d", dest="directory", type=str, required=True)
    return parser.parse_args()

def find_paired_reads(path, ext=".fq.gz"):
    all_reads = sorted([os.path.join(path, f) for f in os.listdir(path) if f.endswith(ext)])
    n_samples = len(all_reads) // 2
    assert n_samples * 2 == len(all_reads)  # Check there is an even number of files.

    for i in range(n_samples):
        yield all_reads[2*i], all_reads[(2*i) + 1]

def batch_paired_files(paired_files, chunk_size):
    batches = []

    active_size = 0
    active = []
    for pair in paired_files:
        size = get_total_size(*pair)

        if active_size + size > chunk_size:
            batches.append(list(active))

            active = list(pair)
            active_size = size
        else:
            active.extend(pair)
            active_size += size

    batches.append(list(active))
    return batches

def get_total_size(*files) -> float:
    total = 0

    for file in files:
        file = os.path.realpath(file)
        byte_size = os.path.getsize(file)
        total += byte_size / 1e9

    return total

if __name__ == "__main__":
    args = parse_args()

    paired_files = find_paired_reads(args.directory)
    batches = batch_paired_files(paired_files, chunk_size=args.chunk_size)

    for i, batch in enumerate(batches):
        with open("batch%d.txt" % i, "w") as f:
            f.write("\n".join(batch))
