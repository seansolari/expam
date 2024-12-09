#!/usr/bin/env python3

# author: sean.solari@monash.edu
# version: 1.0
# date: Mon 30 Sep, 2024

from argparse import ArgumentParser
import os
from typing import List
from expam.database.config import load_database_config, JSONConfig
from expam.tree.tree import Index, Location


def parse_args():
    parser = ArgumentParser(description="Calculate breadth-first enumeration of nodes in a (Newick format) tree.")
    parser.add_argument("-t", "--tree", dest="tree", type=str, help="Path to Newick file")
    parser.add_argument("--db", dest="database", type=str, help="Path to Expam database containing Newick tree")
    parser.add_argument("-o", "--out", dest="out", type=str, help="Path to save output")
    return parser.parse_args()

def breadth_first(tree: Index):
    stack: List[Location] = [tree.pool[1]]
    while stack:
        next_node = stack.pop(0)
        for child in next_node.children:
            stack.append(tree[child])
        yield next_node

if __name__ == "__main__":
    args = parse_args()
    
    # establish tree path
    if args.tree is None and args.database is None:
        print("[error] Must supply either tree (-t, --tree) or Expam database (--db).")
        print("[error] Run --help for more information.")
        exit(1)
    elif args.tree is not None and args.database is not None:
        print("[error] Both tree and database have been supplied. Only supply one.")
        exit(1)
    elif args.tree is not None:
        if not os.path.exists(args.tree):
            print("[error] Could not find input path %s" % args.tree)
            exit(1)
        phylogeny_path = args.tree
    elif args.database is not None:
        if not os.path.exists(args.database):
            print("[error] Could not find input path %s" % args.database)
        database_config = load_database_config(args.database)
        config = JSONConfig(database_config.conf)
        phylogeny_path = config.get_phylogeny_path()

    # load tree
    print("[info] Loading tree from %s" % phylogeny_path)
    _, index = Index.load_newick(phylogeny_path)
    
    # establish output path
    if not phylogeny_path.endswith(".nwk"):
        print("[error] Newick file must end in '.nwk'.")
        exit(1)

    if args.out is None:
        out_file = os.path.basename(phylogeny_path)[:-4] + "-breadthfirst.txt"
    else:
        out_file = args.out

    # breadth-first order
    print("[info] Writing output to %s" % out_file)
    with open(out_file, "w") as f:
        for node in breadth_first(index):
            if node.type == "Branch":
                f.write("p%s\n" % node.name)
            else:
                f.write("%s\n" % node.name)
