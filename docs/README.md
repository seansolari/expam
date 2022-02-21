# expam documentation

All *expam* commands take the following form:
```console
> expam [command] [args...]
```
where `command` is one of several commands detailed below, and `args` are any collection of those outlined in the
**Flags and Options** section. 
Note that some commands require at least a few arguments.

## Commands

### Setting default database.
```console
> expam default_db --db [database name] >> ~/.bash_profile
```
A call to `default_db` will print the required environment variable declaration that *expam* will try to read
at runtime;
```console
> expam default_db --db [database name]
export EXPAM_DEFAULT_DB=[database directory]
```
The >> operator simply redirects this to your bash configuration. This will only be activated in a new session, or 
by running
```console
> source ~/.bash_profile
```
After setting this environment variable, you will no longer need to supply the *--db* flag if you wish to use this
database.
<hr>

### Initialising a database

A new database can be created using the `create` command. In the following, we create a database called
*new* in the current directory.

```console
> expam create -db new
```

#### Build parameters

Database parameters can be set using the `set` command. Parameters to be set for the database build
include
- **-k** or **--kmer**,
- **-n** or **--n_processes**,
- **-p** or **--phylogeny**; Newick format phylogeny.

```console
> expam set -db new -k 31 -n 12 -p ~/my_trees/new.nwk
```

#### Adding sequences

Sequences can be added to the database using the `add` command. This command requires
- **-db** or **--db_name**; the database to which sequences will be added,
- **-d** or **--directory**; the directory to the sequence file, or folder of sequences,
- <u>optionally</u> **--group**; the group to which this sequence will be added (see following note).
- <u>optionally</u> **--first_n**; add first *n* sequences from this directory.

```console
> expam add -db new -d ~/my_sequences/ --first_n 100
```

In the exact sample way, sequences can be removed from the configuration.

```console
> expam remove -db new -d ~/my_sequences/
```

Note that adding and removing sequences ***only impacts future database builds***. Running the 
add and remove commands after a database has been built does not impact the current database.

#### Sequence groups and tree building.

This is a short note on grouping sequences, and the utility of this for building the distance tree (if you are not 
supplying your own). A more in-depth tutorial can be found in [Tutorial 3](./tutorials/3_tree).

Consider the case of building a distance tree containing both viruses and bacteria, with distance estimation
done using *mash*. Trying to estimate genome distance using a single **sketch** and **k-mer** size is unlikely
to be accurate. 

To facilitate tree building, we therefore allow sequences to be added into the database within groups. Say we have
some *bacterial* sequences at `/home/sequences/bacteria/` and *viral* sequences at `/home/sequences/virus/`, 
these can be added to two separate groups;

```console
> expam add -db tutorial -d /home/sequences/virus/ --group virus
...
> expam add -db tutorial -d /home/sequences/bacteria --group bacteria
...
```

Each group can have different values for **k-mer** and **sketch** size. In our example, to accommodate for the
difference in relative genome size between these two groups, we set appropriately scaled parameters;

```console
> expam set -db tutorial --group virus -k 16 -s 1000
...
> expam set -db tutorial --group bacteria -k 21 -s 50000
...
```

If group parameters are not set, *expam* will use the global parameters (those set without the **--group** flag).

*expam* can then run mash to build distance trees for each of these groups individually. We then need to tell
*expam* how to combine these groups, so that a distance tree containing all sequences can be made. Branch length
does not affect the database, only the topological structure. In the case of two groups, their concatenation 
entirely arbitrary,

```console
> cat '({{virus}},{{bacteria}});' > ~/tutorial/phylogeny/tree.nwk
> expam set -db tutorial -p ~/tutorial/phylogeny/tree.nwk
```

Note how we specify groups names using **{{group}}** within the Phylip format.

There are a few options for choosing how *expam* then builds the tree. One option is to feed mash
distances into a local installation of *RapidNJ* to build the tree. Alternatively, *QuickTree* can be used
in place of *RapidNJ* with the appropriate flag;

```console
> expam tree -db tutorial
> expam tree -db tutorial --quicktree
```

Alternatively, a local installation of *mashtree* can be used to sketch sequences and build the tree
in one step.

```console
> expam mashtree -db tutorial
```

A more in-depth walk-through can be found in [Tutorial 3](./tutorials/3_tree).

#### Printing database parameters

A report of the current database configuration can be found through the **print** command.

```console
> expam print -db test
<<< expam configuration file: test >>>

phylogeny       -->     None
k               -->     31
n               -->     4
sketch          -->     None
pile            -->     512

----------------
group name: default
        k               -->     None
        sketch          -->     None
        sequences       -->     9
```

<hr>

### Building the database.

Once all required database parameters (**n_processes**, **phylogeny**, **k-mer size**) have been set,
the database can be build using the **build** command.

```console
> expam build -db tutorial
```

#### Quick build

In the event all database sequences are contained in one folder and one has a distance tree of these
sequences, the database can be initialised and built in the **quickbuild** command.

```console
> expam quickbuild -db my_db -k 31 -n 22 -p ~/trees/my_tree.nwk -d ~/sequences/
```


### Classifying metagenomic reads

Metagenomic reads (in fasta, fastq and compressed *.gz* formats) can be run against the database using
the **run** command. All classification results will be stored in a `results` folder in the database directory.
By default, each *run* command will create a new folder within the results directory named
`.../YYYYMMDD_#/` where *#* corresponds to the id of this run. 

Files containing reads can be run once at a time;

```console
> expam run -db my_db -d ~/reads/reads_1.fna
```

Reads from this run will be located at `~/my_db/results/YYYYMMDD_1`. Reads in some directory can also be
run at the same time,

```console
> expam run -db my_db -d ~/reads/
```

with the results from this run being stored in `~/my_db/results/YYYYMMDD_2`. You can alternatively specify the 
name of the results folder;

```console
> expam run -db my_db ~/reads/ --name my_reads
```

with results stored in `~/my_db/results/my_reads/`. An outline of database runs and results 
can be found in [Tutorial 2](./tutorials/2_classify).

#### Visual depiction of results

By default, running the `to_taxonomy` flag also attempts to overlay the classification results on the tree
used to build the database. The visual features of this tree can be customised.
See [Tutorial 4](./tutorials/4_graphical) for details.

#### Grouping of samples

When running multiple samples, one may wish to combine the counts from a subset of these samples. The `--group` flag
can also be used here, which within the context of a `run` command implies that sample names within this group
will be combined. 

```console
> expam run -db my_db -d ~/reads/ --group sample_1 sample_2 sample_3 --group sample_4 sample_5
```

In this instance, `~/reads/` would contain files that look like `sample_1.fasta` and `sample_4.fa.tar.gz`. 
The samples 1, 2 and 3 will then appear in the results as `sample_1` (the first in your list of groups) and 
samples 4 and 5 will be combined under `sample_4`. 

<hr>

### Converting results to taxonomy

If no additional flags are supplied, classification results will be given with respect to the tree that was used
to make the database. These results can be found in the `phy` subdirectory of your results folder.

#### Downloading taxonomy data

You may also be interested in a taxonomic overview of classification results. To enable this, *expam* first
needs to download the associated taxonomy for all sequences in the database. Taxonomy data by default is found
by querying the NCBI Accession ID from the header of each sequence through the NCBI Entrez platform. Alternatively, 
you can also explicitly provide a taxonomic id for any given sequence by placing `taxid|[taxid]` at the header
line of the sequence. For example,

```text
>NZ_CP063205.1|taxid|35746
```

Taxonomy data is downloaded using the `download_taxonomy` command. ***This can only be run after a database has
been built.*** 

```console
> expam download_taxonomy -db my_db
```

*expam* will only download the required data for the sequences in your database, not the entire NCBI taxonomy. This 
has the advantage of saving disk space, and is a convenient command to keep the taxonomic information up-to-date.

#### Taxonomic summary of results

After running the `download_taxonomy` command, you now have access to both the `--taxonomy` flag which can be 
used with the `run` command, as well as the `to_taxonomy` command. 

Supplying the `--taxonomy` flag for a classification `run` will immediately provide both phylogenetic and taxonomic
summary of the results.

```console
> expam run -db my_db ~/reads/ --name my_run_name --taxonomy
```

Alternatively, phylogenetic results can be converted to taxonomic results after-the-fact.

```console
> expam to_taxonomy -db my_db --name old_run_name
```

<hr>


## Flags and Options

Commands preceded by [ ! ] are more important to know than the others. Example usages can 
be found in the above section, and also in the tutorial pages.

<hr>

#### [ ! ] -db, --db_name
Directory of database on which command will be supplied.

<hr>

#### [ ! ] -k, --kmer
Length of k for k-mers used in the analysis.

<hr>

#### [ ! ] -n, --n_processes
Number of Python processes to be used for the analysis.

<hr>

#### -s, --sketch
Sketch size to be used for tree building. See URL for examples.

<hr>

#### [ ! ] -p, --phylogeny
Directory of phylogeny to be used for building the database.

<hr>

#### [ ! ] -d, --directory
Used to supply the directory for some files that are required for a command.

<hr>

#### --sourmash
Use sourmash for distance estimation and genome sketching. See URL for examples.

<hr>

#### --rapidnj
Use a local installation of RapidNJ () for Neighbour-Joining the distance matrix.

<hr>

#### --quicktree
Use a local installation of QuickTree () for Neighbour-Joining the distance matrix.

<hr>

#### -l, --length
Length of reads to be sampled. Useful for unit testing.

<hr>

#### -o, --out
Used to supply the directory where the results of some analysis will be saved.

<hr>

#### -y, --pile
Number of genomes to pile at a time when building the database. Must be of the form
2^K for some integer K.

<hr>

#### -e, --error_rate
Generate errors in uniformly sampled reads (error_rate ~ reads with errors / reads)

<hr>

#### --plot
Plot the timing graph after database building completion.

<hr>

#### --first
Supply a cutoff for the number of sequences to be added from some directory. Useful for unit testing.

<hr>

#### [ ! ] --cutoff
Ignore organisms in the results with less that `cutoff` reads. Useful to ignore reads that have
been affected by read error.

<hr>

#### [ ! ] --cpm
Cutoff value specified in counts per million.

<hr>

#### [ ! ] --taxonomy
Convert phylogenetic results to taxonomic results. This is to be used with the `run` command.

<hr>

#### --phyla
Colour the four main Gut phyla on the output phylotree.

<hr>

#### --keep_zeros
Keep nodes in the output where no reads have been assigned.

<hr>

#### --ignore_names
Don't plot genome names on the outer ring of the phylotree.

<hr>

#### --group
*Use depends on context.*

With commands `run` or `to_taxonomy`, this is a space separated list, potentially starting with
a colour in hexadecimal form, of sample names whose results are to be pooled together in the phylotree.
See [Tutorial 4](./tutorials/4_graphical) for examples.

With database initialisation commands `add`, `remove` or `set`, this specifies which group of sequences
is to have their parameters modified. See [Tutorial 1](./tutorials/1_basic) for examples.

<hr>

#### --colour_list
List of colours, which will be used in the order provided, to be used when choosing the colours for plotting
samples/groups in the phylotree.

<hr>

#### --name
Name of output results folder.

<hr>

#### --circle_scale
Scale the size of circles plotted atop the phylotree which represent reads classified as splits.

