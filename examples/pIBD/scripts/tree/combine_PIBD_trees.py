#!/usr/bin/env python3
import re

def read_newick(path):
    with open(path, 'r') as f:
        data = f.read().strip().rstrip(";")

    data = data.replace(".fna", "")
    return data

def complete_template(template: str):
    subtrees = re.findall(r"{\S+}", template)
    _format = lambda s: s[1:-1]
    newick_data = tuple(read_newick(_format(file)) for file in subtrees)
    
    final_template = template
    for subtree in subtrees:
        final_template = final_template.replace(subtree, "%s")

    return final_template % newick_data

if __name__ == "__main__":
    template = "(({bacteria.nwk}, {archaea.nwk}), {virus.nwk});"
    print(complete_template(template))
