Database Build
==============

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

.. warning::

    When using custom reference sequences or assembled genomes, ensure as 
    much as possible that these genomes are uncontaminated - contaminated
    genomes compromise distance estimation for building the reference
    tree and subsequently may impact classification performance.

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
* We'll first ensure that :code:`sourmash` is installed, before running the :code:`tree` command to build the tree.

.. code-block:: console

    $ python3 -m pip install sourmash
    $ expam tree -db test --sourmash
    ...

* :code:`print` and match with the following output:

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


Build the database
------------------

* Now that we have a distance-tree for our added reference sequences, we can build the database.

.. code-block:: console

    $ expam build -db test
    Clearing old log files...
    Importing phylogeny...
    * Initialising node pool...
    * Checking for polytomies...
        Polytomy (degree=3) detected! Resolving...
    * Finalising index...
    Creating LCA matrix...
    Extracting sequences from /Users/ssol0002/Documents/Projects/pam/test/data/sequences/GCF_000008725.1_ASM872v1_genomic.fna.gz...
    Extracting sequences from /Users/ssol0002/Documents/Projects/pam/test/data/sequences/GCF_000007765.2_ASM776v2_genomic.fna.gz...
    Extracting sequences from /Users/ssol0002/Documents/Projects/pam/test/data/sequences/GCF_000005845.2_ASM584v2_genomic.fna.gz...
    Extracting sequences from /Users/ssol0002/Documents/Projects/pam/test/data/sequences/GCF_000006925.2_ASM692v2_genomic.fna.gz...
    Extracting sequences from /Users/ssol0002/Documents/Projects/pam/test/data/sequences/GCF_000006945.2_ASM694v2_genomic.fna.gz...
    Extracting sequences from /Users/ssol0002/Documents/Projects/pam/test/data/sequences/GCF_000006765.1_ASM676v1_genomic.fna.gz...
    expam: 42.359643852s

    PID - 65856 dying...