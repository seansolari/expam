## Building a reference phylogeny

We have a collection of reference genomes from Archaea, Bacteria and Viruses. We will build a mash tree for each of these domains, and then concatenate these trees to form the reference phylogeny for *expam*.

### Sketching the reference genomes.

The reference genomes are located in `bacteria/`, `archaea/` and `virus/` in FASTA format. Sketch the sequences with `mash`.

```bash
mash sketch -k 21 -s 50000 -p 32 -o bacteria.k21.s50000 bacteria/*.fna
mash sketch -k 21 -s 50000 -p 10 -o archaea.k21.s50000 archaea/*.fna
mash sketch -k 16 -s 5000 -p 10 -o virus.k16.s5000 virus/*.fna
```

### Construct distance matrices and format the output.

Compute pairwise distances

```bash
mash dist -p 32 -t *.msh *.msh > distances.tab
```

Format the results for `rapidnj`.

```bash
python3 scripts/tree/format_distance_matrix.py distances.tab > formatted_distances.tab
```

### Run Neighbour-Joining

Run `rapidnj`...

```bash
rapidnj formatted_distances.tab -i pd -o t -c 32 > tree.nwk
```

...and format the output Newick file to clean genome names.

```bash
python3 scripts/tree/format_rapidnj_tree.py tree.nwk > new_tree.nwk
```

### Combine trees from each domain

Combine trees using the Newick template `(({bacteria}, {archaea}), {virus});`

```bash
python3 scripts/tree/combine_PIBD_trees.py > pibd.nwk
```

### Create the *expam* database

Build an expam database from these reference sequences.

```bash
expam create -db pibd
expam set -db pibd -k 31 -n 32 -p pibd.nwk
expam add -db pibd -d archaea/
expam add -db pibd -d bacteria/
expam add -db pibd -d virus/
expam build -db pibd
```

**Next step...**

[Batching samples, removing human reads and classifying](./1_classification.md)
