Quickstart
==========

Creating a database
-------------------
* We create a database called :code:`test`.

.. code-block:: console

    $ expam create -db test
    Successfully created directories /Users/seansolari/Documents/Databases/test/phylogeny!
    Fresh database configuration generated!
    Logs path /Users/seansolari/Documents/Databases/test/logs created!
    Made path /Users/seansolari/Documents/Databases/test/database.
    Made path /Users/seansolari/Documents/Databases/test/results.

* This database will be located in the current directory, in a folder called :code:`test`.


Add reference sequences
^^^^^^^^^^^^^^^^^^^^^^^
We now add sequences to our new database.

* We have supplied sequences with the source:
  
    * Reference genomes: :code:`../expam/test/data/sequence/`,
    * Reads (for later): :code:`../expam/test/data/reads/`.

* Add these sequences to our database :code:`test`.

.. code-block:: console

    $ expam add -db test -d ~/Documents/expam/test/data/sequences
    Added 6 files from /Users/seansolari/Documents/expam/test/data/sequences/.

.. note::

    **expam** can handle both compressed and uncompressed sequences files and
    automatically detects how it should open most regular sequence files.


Specify build parameters
^^^^^^^^^^^^^^^^^^^^^^^^
* These are bacterial genomes, so :code:`k=31` should work fine.
* We'll use 4 processes to speed up the build.

.. code-block:: console

    $ expam set -db test -k 31 -n 4


Build a tree for the reference database
---------------------------------------

* Expam needs a tree with the phylogenetic relationship between your sequences.
* We'll make a tree out of :code:`mash` distances.

* :code:`sourmash` is a portable Python version of :code:`mash` (you can install it with `pip <https://pypi.org/project/sourmash/>`_).
* We'll use a sketch size of :code:`s=1000` to represent each genome, and then compare them.

.. code-block:: console

    $ expam set -db test -s 1000

* We'll use :code:`RapidNJ` to make a tree from the :code:`sourmash` distances (see :doc:`here <../dependencies>` to install).
* Run the :code:`tree` command to build the tree.

.. code-block:: console

    $ expam tree -db test --sourmash
    ...

* :code:`print` and match with my output.

.. code-block:: console

    $ expam print -db test
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


Running classifications
-----------------------

* Recall the reads we supplied:
  
    * :code:`../expam/test/data/reads/`
  
* We use the :code:`run` command to classify reads.
* These are paired reads, but for now we'll treat them as separate.
* By default, run results are stored in the :code:`results` database folder,
  
    * here :code:`test/results`.
    * This can be redirected using :code:`--out`.
  
* We can supply a :code:`--name` to label these results.
  
    * We'll call this first run :code:`unpaired`.

.. code-block:: console
    
    $ expam run -db test -d /Users/seansolari/Documents/expam/test/data/reads/ --name unpaired_test
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

* The results can be found at :code:`test/results/unpaired_test/`.

.. code-block:: console

    $ ls test/results/unpaired_test
    phy     phylotree.pdf
    $ ls test/results/unpaired_test/phy/
    GCF_000005845.2_ASM584v2_genomic.gz_1.csv       classified_counts.csv                           splits_counts.csv
    GCF_000005845.2_ASM584v2_genomic.gz_2.csv       raw


classified_counts.csv
^^^^^^^^^^^^^^^^^^^^^

* This file gives an overview of the confident classification results.

.. code-block:: console

    $ cat test/results/unpaired_test/phy/classified_counts.csv
            GCF_000005845.2_ASM584v2_genomic.gz_2   GCF_000005845.2_ASM584v2_genomic.gz_1
    unclassified    0       0
    p1      3       3
    p2      232     232
    GCF_000005845.2_ASM584v2_genomic        765     765


GCF_000005845.2_ASM584v2_genomic.gz_1.csv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Each input sample is provided with its own summary file, with the following (tab-delimited) columns:
  
    * Phylogenetic node,
    * % confident cumulative,
    * classified cumulative (at or below this node),
    * raw classification count (at this node),
    * % split cumulative,
    * split cumulative,
    * raw split count.

.. code-block:: console

    $ cat test/results/unpaired_test/phy/GCF_000005845.2_ASM584v2_genomic.gz_1.csv
    unclassified    0.000000%       0       0                       
    p1      100.000000%     1000    3       0.000000%       0       0
    p2      99.700000%      997     232     0.000000%       0       0
    GCF_000005845.2_ASM584v2_genomic        76.500000%      765     765     0.000000%       0       0


raw/GCF_000005845.2_ASM584v2_genomic.gz_1.csv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Read-wise output for all input samples can be found in the :code:`raw` subdirectory.

.. code-block:: console

    $ head test/results/unpaired_test/phy/raw/GCF_000005845.2_ASM584v2_genomic.gz_1.csv
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


Running paired data
-------------------

* To run paired-end data, supply the :code:`--paired` flag.

.. note::

    **expam** automatically searches the input files and matches those most similar file names as paired.

* We'll call this separate run :code:`paired_test`.

.. code-block:: console

    $ expam run -db test -d /Users/seansolari/Documents/expam/test/data/reads/ --name paired_test --paired
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


Taxonomic results
-----------------

