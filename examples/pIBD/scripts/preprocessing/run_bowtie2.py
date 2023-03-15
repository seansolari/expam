#!/usr/bin/env python3
from argparse import ArgumentParser
import os
import subprocess
import sys

def read_batch(file):
    with open(file, 'r') as f:
        batch_members = [line.strip() for line in f]

    batch_members = sorted(batch_members)
    n_samples = len(batch_members) // 2
    assert n_samples * 2 == len(batch_members)  # Check there is an even number of files.

    for i in range(n_samples):
        yield batch_members[2*i], batch_members[(2*i) + 1]

def string_intersect(s1: str, s2: str) -> str:
    for i, (a, b) in enumerate(zip(s1, s2)):
        if a != b:
            return s1[:i]
    else:
        m = min(len(s1), len(s2))
        return s1[:m]

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--batch", dest="batch_file", type=str, required=True)
    parser.add_argument("-o", dest="out_path", type=str, required=True)
    return parser.parse_args()

def run_bowtie(forward, reverse, out_dir, threads=12):
    forward = os.path.abspath(forward)
    reverse = os.path.abspath(reverse)
    out_dir = os.path.abspath(out_dir)

    out_file = string_intersect(os.path.basename(forward), os.path.basename(reverse)).rstrip(".") + ".unmapped.fq.gz"
    out_path = os.path.join(out_dir, out_file)

    bowtie_command = [
        "bowtie2",
        "-p", str(threads),
        "-x", "GRCh38_noalt_as",
        "-1", forward, "-2", reverse,
        "--un-conc-gz", out_path
    ]

    p = subprocess.run(bowtie_command, stdout=sys.stdout)
    if p.stderr:
        print("Error encountered:")
        for line in p.stderr:
            print("\t" + line.strip())

        sys.exit(p.returncode)

    print("Unmapped reads written to %s.\n" % out_path)

if __name__ == "__main__":
    args = parse_args()

    for forward, reverse in read_batch(args.batch_file):
        print("Running bowtie2 on\n\tF=%s\n\tR=%s" % (forward, reverse))
        run_bowtie(forward, reverse, args.out_path)
