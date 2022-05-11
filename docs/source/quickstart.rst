Quickstart Tutorial
===================

We will be using a pre-built database to classify some metagenomic reads and obtain both phylogenetic and taxonomic output.

Get the database
----------------

* Download one of the following compressed databases (see :doc:`Tutorial 1 <tutorials/overview>` for instructions on building a database).

.. container:: wideflexcontainer

   .. container:: colone

      .. raw:: html

         <p>
            <a class="reference internal" href="https://drive.google.com/file/d/1K1sVA4LGgmGBVg_0GeUppVa_xcfWxxGL/view?usp=sharing" target="_blank">Test Database (110.7 Mb)</a>
         </p>

   .. container:: coltwo

      .. raw:: html

         <p>
            <a class="reference internal" href="https://figshare.com/s/3475c3a9aa926a40c722" target="_blank">expam RefSeq (122.35 Gb)</a>
         </p>

* Unzip the compressed file. For instance, if you downloaded :code:`Test Database` (:code:`test.tar.gz`) into your home directory :code:`~`, run

.. code-block:: console

    $ tar -xvzf test.tar.gz

The database will now be located at :code:`~/test`, which is the directory you should pass to any :code:`-db` flag as input to **expam**.

* Download these metagenomic reads (simulated reads from a genome used to build the :code:`test` database).

.. container:: wideflexcontainer

   .. container:: colone

      .. raw:: html

         <p>
            <a class="reference internal" href="https://drive.google.com/file/d/1hTOndUelxf1cEEW8EIRYZrdxaHKez9qz/view?usp=sharing" target="_blank">Fasta reads (432 Kb)</a>
         </p>

* Unzip these metagenomic reads into a new folder, which we will call :code:`reads`. Assuming you have downloaded and moved the above reads into your home directory, run

.. code-block::

    $ mkdir reads
    $ mv reads.tar.gz reads
    $ cd reads
    $ tar -xvzf reads.tar.gz

There should be two files:

    * :code:`~/reads/GCF_000005845.2_ASM584v2_genomic.fna.gz_1.fa`,
    * :code:`~/reads/GCF_000005845.2_ASM584v2_genomic.fna.gz_2.fa`.

These are paired read, fasta files.


Phylogenetic classification
---------------------------

* We will now classify these reads using the database you downloaded. We will save the results to an output folder located at :code:`~/my_run/`.

* Run the :code:`expam classify` command as follows (replacing :code:`~/test` with where you decompressed the database from Step 1 if necessary):

.. code-block:: console

    $ expam classify -db ~/test -d ~/reads --paired --out ~/my_run
    Clearing old log files...
    Loading the map and phylogeny.

    Preparing shared memory allocations...
    Loading database keys...
    Loading database values...
    * Initialising node pool...
    * Checking for polytomies...
        Polytomy (degree=3) detected! Resolving...
    * Finalising index...
    Loading reads from /Users/ssol0002/Documents/Projects/pam/test/data/reads/GCF_000005845.2_ASM584v2_genomic.fna.gz_2.fa, /Users/ssol0002/Documents/Projects/pam/test/data/reads/GCF_000005845.2_ASM584v2_genomic.fna.gz_1.fa...
    Could not import ete3 plotting modules! Error raised:
    Traceback (most recent call last):
    File "/Users/ssol0002/Documents/Projects/pam/src/expam/tree/tree.py", line 622, in draw_tree
        import ete3.coretype.tree
    ModuleNotFoundError: No module named 'ete3'

    Skipping plotting...
    Could not import ete3 plotting modules! Error raised:
    Traceback (most recent call last):
    File "/Users/ssol0002/Documents/Projects/pam/src/expam/tree/tree.py", line 622, in draw_tree
        import ete3.coretype.tree
    ModuleNotFoundError: No module named 'ete3'

    Skipping plotting...

.. note::

    Note that **expam** tried to plot the results on a phylotree, but since we did not have the ete3 module installed,
    it simply skipped plotting the results. This is the expected behaviour to let you know **expam** was not able
    to produce a graphical picture for your results.


* The phylogenetic classifications will be located at :code:`~/my_run/phy`, and will contain four files:
    * :code:`~/my_run/GCF_000005845.2_ASM584v2_genomic.gz_1.csv` - sample summary file,

    .. code-block::

        unclassified    0.000000%       0       0                       
        p1      100.000000%     1000    3       0.000000%       0       0
        p2      99.700000%      997     232     0.000000%       0       0
        GCF_000005845.2_ASM584v2_genomic        76.500000%      765     765     0.000000%       0       0

    * :code:`~/my_run/classified.csv` - classified summary file,

    .. code-block::

                GCF_000005845.2_ASM584v2_genomic.gz_1
        unclassified    0
        p1      3
        p2      232
        GCF_000005845.2_ASM584v2_genomic        765

    * :code:`~/my_run/split.csv` - split summary file,

    .. code-block:: 

                GCF_000005845.2_ASM584v2_genomic.gz_1
        p1      0
        p2      0
        GCF_000005845.2_ASM584v2_genomic        0

    * :code:`~/my_run/raw` - raw read-wise classifications. There will be a single raw read-wise output file, :code:`~/my_run/raw/GCF_000005845.2_ASM584v2_genomic.gz_1.csv`.

    .. code-block::

        C       R4825323246286034638    p2      302     p2:240
        C       R4280015672552393909    p10     302     p10:240
        C       R5925738157954038177    p10     302     p1:5 p10:16 p2:198 p10:16 p1:5
        C       R3237657389899545456    p10     302     p2:85 p10:31 p2:8 p10:31 p2:85
        C       R6111671585932593081    p10     302     p2:36 p10:37 p2:3 p10:88 p2:3 p10:37 p2:36
        C       R4574482278193488645    p10     302     p10:29 p2:14 p10:31 p2:2 p10:88 p2:2 p10:31 p2:14 p10:29
        C       R8975058804953044791    p10     302     p10:21 p2:59 p10:80 p2:59 p10:21
        C       R6052336354009855322    p10     302     p2:53 p10:31 p2:72 p10:31 p2:53

The sample summary file is a tab-separated document where the first element of each row is a phylogenetic node/clade, and the corresponding values are contain details of the raw and cumulative classifications and splits at this particular node.

The classified summary file is a tab-separated matrix where each row is a phylogenetic clade, each column is an input sample, and the cell value is the raw counts at this clade. The split summary file is an analogous file that contains the raw split count at any given clade. These two files are formatted such that they will always have the same column and row indices, and in the same order.

The raw read-wise output is a sub-directory containing one output file for each input sample, the kraken-formatted read-wise output.

A more comprehensive overview is given :doc:`this tutorial <tutorials/classify>`.


Convert to taxonomy
-------------------

* First run :code:`expam download_taxonomy` download the taxonomy for all sequences in the database. This will require an internet connection.

.. code-block:: console

    $ expam download_taxonomy -db ~/test
    Posting 6 UIDs to NCBI Entrez nuccore.
    Received 6 response(s) for ESummary TaxID request!
    Posting 6 UIDs to NCBI Entrez taxonomy.
    Received 6 response(s) for EFetch Taxon request!
    Taxonomic lineages written to ~/test/phylogeny/taxid_lineage.csv!
    Taxonomic ranks written to ~/test/phylogeny/taxa_rank.csv!

* We saved our previous classification results to :code:`~/my_run`. This is the directory we pass to :code:`expam to_taxonomy` to convert phylogenetic classifications to taxonomy.

.. code-block:: console

    $ expam to_taxonomy -db test --out ~/my_run
    * Initialising node pool...
    * Checking for polytomies...
        Polytomy (degree=3) detected! Resolving...
    * Finalising index...


* There will now be taxonomic output files located in :code:`~/my_run/tax/`, analogous to each of the files present in the phylogenetic output, with the exception of :code:`classified.tsv` and :code:`split.tsv` - only the sample summaries and raw read-wise output are converted.

    * :code:`~/my_run/tax/GCF_000005845.2_ASM584v2_genomic.gz_1.csv` - taxonomic summary file

    .. code-block::

                c_perc  c_cumul c_count s_perc  s_cumul s_count rank    scientific name
        unclassified    0.0%    0       0       0%      0       0       0       0
        1       100.0%  1000    0       0%      0       0       root    
        131567  100.0%  1000    0       0%      0       0       top     cellular organisms
        2       100.0%  1000    235     0%      0       0       superkingdom    cellular organisms Bacteria
        1224    76.5%   765     0       0%      0       0       phylum  cellular organisms Bacteria Proteobacteria

    * :code:`~/my_run/tax/raw/GCF_000005845.2_ASM584v2_genomic.gz_1.csv` - taxonomic read-wise output. The second column is the read header, the third is the assigned taxid, and the fourth is the length of the read. Observe length of 300 for paired-end 150bp reads.

    .. code-block::

        C       R4825323246286034638    2       302
        C       R4280015672552393909    511145  302
        C       R5925738157954038177    511145  302
        C       R3237657389899545456    511145  302
        C       R6111671585932593081    511145  302
        C       R4574482278193488645    511145  302
        C       R8975058804953044791    511145  302
        C       R6052336354009855322    511145  302
        C       R4978825024774141837    2       302
        C       R7016203356160788326    511145  302

The complete comprehensive overview is given :doc:`this tutorial <tutorials/classify>`.

