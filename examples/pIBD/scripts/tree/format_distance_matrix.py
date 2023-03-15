#!/usr/bin/env python3
import os
import sys

def format_matrix(file):
    with open(file, 'r') as f:
        names = f.readline().strip().split("\t")
        print("%d" % len(names[1:]))

        for line in f:
            parts = line.strip().split("\t")
            parts[0] = os.path.basename(parts[0])
            print("\t".join(parts))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Invalid arguments!")
        sys.exit(1)

    file = sys.argv[1]
    format_matrix(file)

