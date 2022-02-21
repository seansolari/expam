# 1. A 'Hello World' example.

## Creating the database.
We will create a new database, called `test`.
```bash
> expam create -db test
Successfully created directories /Users/seansolari/Documents/Databases/test/phylogeny!
Fresh database configuration generated!
Logs path /Users/seansolari/Documents/Databases/test/logs created!
Made path /Users/seansolari/Documents/Databases/test/database.
Made path /Users/seansolari/Documents/Databases/test/results.
```

### Add Sequences
Within the source files for expam, there are sequences located at `.../expam/test/data/sequences/`. 
These sequences are from the RefSeq collection, and are all compressed `.fna.gz` files. `expam` can handle
both compressed and uncompressed sequences files. We will use these sequences in our database build,
so we should make the `test` database aware of them.

```bash
> expam add -db test -d ~/Documents/expam/test/data/sequences/
Added 6 files from /Users/seansolari/Documents/expam/test/data/sequences/.
```

By specifying a folder, expam will add all sequence files immediately contained in that directory.

### Specify build parameters.
To build this database, we will use a k-mer length of `31`, and a total of `4` processes, among which
the building work will be distributed. The more processes you use, the shorter the build time. 
```bash
> expam set -db test -k 31 -n 4 
```

### Build A Tree
To build a database, expam requires a tree that specifies the relationship between your sequences. 
To be most effective, sequences that are closer together in the tree should have similar sequences. For
this reason, we build a tree based on mash distances. 

We will be using a Python variant of mash, `sourmash`, to build our tree. See the <b>4_tree</b> tutorial
for a more in-depth guide to building a tree. Here we will use a sketch size of `1000`, which is likely 
too small for bacterial sequences, but is computationally light.

```bash
> expam set -db test -s 1000
```

This tutorial uses RapidNJ for building the tree after using sourmash to sketch sequence distances. If you
do not have an installation of RapidNJ, a convenient method to install it would be to create a local
conda environment, and run the following.

```bash
> conda install rapidnj
```

We can now run the `tree` command to build the tree.
```bash
> expam tree -db test --sourmash
...
> expam print -db test
<<< expam configuration file: test >>>

phylogeny       -->     /Users/seansolari/Documents/Databases/test/phylogeny/tree/test.nwk
k               -->     31
n               -->     4
sketch          -->     1000
pile            -->     None

----------------
group name: default
        k               -->     None
        sketch          -->     None
        sequences       -->     6

```

You should see a similar output after running the `print` command.

### Build the database.
The database can now be build, using the `build` command.
```bash
> expam build -db test

Clearing old log files...
Importing phylogeny...
* Initialising node pool...
* Checking for polytomies...
        Polytomy (degree=3) detected! Resolving...
* Finalising index...
Creating LCA matrix...
Extracting sequences from /Users/seansolari/Documents/expam/test/data/sequences/GCF_000008725.1_ASM872v1_genomic.fna.gz...
Extracting sequences from /Users/seansolari/Documents/expam/test/data/sequences/GCF_000007765.2_ASM776v2_genomic.fna.gz...
Extracting sequences from /Users/seansolari/Documents/expam/test/data/sequences/GCF_000005845.2_ASM584v2_genomic.fna.gz...
Extracting sequences from /Users/seansolari/Documents/expam/test/data/sequences/GCF_000006925.2_ASM692v2_genomic.fna.gz...
Extracting sequences from /Users/seansolari/Documents/expam/test/data/sequences/GCF_000006945.2_ASM694v2_genomic.fna.gz...
Extracting sequences from /Users/seansolari/Documents/expam/test/data/sequences/GCF_000006765.1_ASM676v1_genomic.fna.gz...
expam: 25.799587115999998s

PID - 92722 dying...
```

## Running classifications.

### Unpaired data
Alongside the collection of reference sequences, there are some simulated reads that can be run against the database.
These have been sampled from one of the reference sequences. These reads are located at `.../expam/test/data/reads/`.

Use the `run` command to run reads against the database. These are actually paired reads, but we will treat them 
separately for the time being. By default, run results are stored in the `results` subdirectory of the database. You
can assign a name to each run, or one will be automatically assigned for you.

We will call this first run `unpaired_test`.

```bash
> expam run -db test -d /Users/seansolari/Documents/expam/test/data/reads/ --name unpaired_test

Clearing old log files...
Results directory created at /Users/seansolari/Documents/Databases/test/results/unpaired_test.
Loading the map and phylogeny.

Preparing shared memory allocations...
Loading database keys...
Loading database values...
* Initialising node pool...
* Checking for polytomies...
        Polytomy (degree=3) detected! Resolving...
* Finalising index...
Loading reads from /Users/seansolari/Documents/expam/test/data/reads/GCF_000005845.2_ASM584v2_genomic.fna.gz_1.fa...
Loading reads from /Users/seansolari/Documents/expam/test/data/reads/GCF_000005845.2_ASM584v2_genomic.fna.gz_2.fa...
Phylogenetic tree written to /Users/seansolari/Documents/Databases/test/results/unpaired_test/phylotree.pdf!
```

The results can be found at `.../test/results/unpaired_test/`.
```bash
> ls test/results/unpaired_test/
phy             phylotree.pdf

> ls test/results/unpaired_test/phy/
GCF_000005845.2_ASM584v2_genomic.gz_1.csv       classified_counts.csv                           splits_counts.csv
GCF_000005845.2_ASM584v2_genomic.gz_2.csv       raw
```

`classified_counts.csv` gives an overview of the confident classification results.
```bash
> cat test/results/unpaired_test/phy/classified_counts.csv
        GCF_000005845.2_ASM584v2_genomic.gz_2   GCF_000005845.2_ASM584v2_genomic.gz_1
unclassified    0       0
p1      3       3
p2      232     232
GCF_000005845.2_ASM584v2_genomic        765     765
```

`GCF_000005845.2_ASM584v2_genomic.gz_1.csv` provides an overview of the classifications for this particular sample.
```bash
> cat test/results/unpaired_test/phy/GCF_000005845.2_ASM584v2_genomic.gz_1.csv
unclassified    0.000000%       0       0                       
p1      100.000000%     1000    3       0.000000%       0       0
p2      99.700000%      997     232     0.000000%       0       0
GCF_000005845.2_ASM584v2_genomic        76.500000%      765     765     0.000000%       0       0
```

The classifications for each read can be found in the `raw` subdirectory.
```bash
> head test/results/unpaired_test/phy/raw/GCF_000005845.2_ASM584v2_genomic.gz_1.csv
C       R4825323246286034638    p2      151     p2:120
C       R4280015672552393909    p10     151     p10:120
C       R5925738157954038177    p10     151     p2:99 p10:16 p1:5
C       R3237657389899545456    p10     151     p2:4 p10:31 p2:85
C       R6111671585932593081    p10     151     p10:44 p2:3 p10:37 p2:36
C       R4574482278193488645    p10     151     p10:44 p2:2 p10:31 p2:14 p10:29
C       R8975058804953044791    p10     151     p10:40 p2:59 p10:21
C       R6052336354009855322    p10     151     p2:36 p10:31 p2:53
C       R4978825024774141837    p2      151     p2:36 p1:34 p2:31 p1:17 p2:2
C       R7016203356160788326    p10     151     p2:4 p10:52 p2:64
```

### Paired data
To run paired-end data, simply supply the `--paired` flag. Expam will automatically try to find all pairs of data. 
We will call this run `paired_test`.

```bash
> expam run -db test -d /Users/seansolari/Documents/expam/test/data/reads/ --name paired_test --paired

Clearing old log files...
Results directory created at /Users/seansolari/Documents/Databases/test/results/paired_test.
Loading the map and phylogeny.

Preparing shared memory allocations...
Loading database keys...
Loading database values...
* Initialising node pool...
* Checking for polytomies...
        Polytomy (degree=3) detected! Resolving...
* Finalising index...
Loading reads from /Users/seansolari/Documents/expam/test/data/reads/GCF_000005845.2_ASM584v2_genomic.fna.gz_2.fa, /Users/seansolari/Documents/expam/test/data/reads/GCF_000005845.2_ASM584v2_genomic.fna.gz_1.fa...
Phylogenetic tree written to /Users/seansolari/Documents/Databases/test/results/paired_test/phylotree.pdf!
```

### Taxonomic results.

Classifications above are given with respect to the phylogeny/tree that was used to create the database. 
There is also an option to give a taxonomic summary of the results.

To do this, we first need to download the associated metadata for all sequences from the database. Expam
automatically collects this data from the NCBI Entrez servers through the `download_taxonomy` command. This can
only be run after the database has been built.

Note that it only downloads the required metadata, not the entire NCBI taxonomy.
This was done to save space.

```bash
> expam download_taxonomy -db test

Posting 6 UIDs to NCBI Entrez nuccore.
Received 6 response(s) for ESummary TaxID request!
Posting 6 UIDs to NCBI Entrez taxonomy.
Received 6 response(s) for EFetch Taxon request!
Taxonomic lineages written to /Users/seansolari/Documents/Databases/test/phylogeny/taxid_lineage.csv!
Taxonomic ranks written to /Users/seansolari/Documents/Databases/test/phylogeny/taxa_rank.csv!
```

We will now convert the previous `paired_test` run results to a taxonomic format. This is done using the `to_taxonomy`
command.

```bash
> expam to_taxonomy -db test --name paired_test

* Initialising node pool...
* Checking for polytomies...
        Polytomy (degree=3) detected! Resolving...
* Finalising index...
Phylogenetic tree written to /Users/seansolari/Documents/Databases/test/results/paired_test/phylotree.pdf!
```

Note now the additional `tax` subdirectory in the results folder, containing the taxonomic summaries.
```bash
> ls test/results/paired_test/tax/
GCF_000005845.2_ASM584v2_genomic.gz_2.csv       raw

> head test/results/paired_test/tax/GCF_000005845.2_ASM584v2_genomic.gz_2.csv
        c_perc  c_cumul c_count s_perc  s_cumul s_count rank    scientific name
unclassified    0.0%    0       0       0%      0       0       0       0
1       100.0%  1000    0       0%      0       0       root    
131567  100.0%  1000    0       0%      0       0       top     cellular organisms
2       100.0%  1000    235     0%      0       0       superkingdom    cellular organisms Bacteria
1224    76.5%   765     0       0%      0       0       phylum  cellular organisms Bacteria Proteobacteria
1236    76.5%   765     0       0%      0       0       class   cellular organisms Bacteria Proteobacteria Gammaproteobacteria
91347   76.5%   765     0       0%      0       0       order   cellular organisms Bacteria Proteobacteria Gammaproteobacteria Enterobacterales
543     76.5%   765     0       0%      0       0       family  cellular organisms Bacteria Proteobacteria Gammaproteobacteria Enterobacterales Enterobacteriaceae
561     76.5%   765     0       0%      0       0       genus   cellular organisms Bacteria Proteobacteria Gammaproteobacteria Enterobacterales Enterobacteriaceae Escherichia
```

All the summary files mentioned above for tree-based results have a counterpart taxonomic result file.
