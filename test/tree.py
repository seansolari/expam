#!/usr/bin/env python3
from expam.tree.tree import Index


def test_tree_reduction():
    base_tree = "(B:1.0,(D:1.0,(F:1.0,G:1.0)E:1.0)C:1.0)A;"
    _, tree = Index.from_newick(base_tree, keep_names=True)
    print(tree.to_newick())
    tree.reduce(["C", "G", "A"])
    print(tree.to_newick())


if __name__ == "__main__":
    test_tree_reduction()
