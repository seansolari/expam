#!/usr/bin/env python3
import sys

def format_tree_string(nwk):
    seq_types = [".fna", ".faa"]
    comp_types = [".tar.gz", ".gz"]

    nwk = nwk.replace("'", "")

    for seq_type in seq_types:
        for comp_type in comp_types:
            ext = seq_type + comp_type + ":"
            nwk = nwk.replace(ext, ":")

    return nwk

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Invalid arguments!")
        sys.exit(1)

    file = sys.argv[1]
    with open(file, 'r') as f:
        nwk_string = f.read().strip()

    print(format_tree_string(nwk_string))

