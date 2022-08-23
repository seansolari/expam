Quickstart Tutorial
===================

We will be using a pre-built database to classify some simulated metagenomic reads and obtain both phylogenetic and taxonomic output.

Get the database
----------------

* Download one of the following compressed databases (see :doc:`Tutorial 1 <tutorials/build>` for instructions on building a database).

.. container:: wideflexcontainer

   .. container:: colone

      .. raw:: html

         <p>
            <a class="reference internal" href="https://drive.google.com/file/d/1KRrEvG5Sr28wvkEFW5CXKwcYqLHAu6f7/view?usp=sharing" target="_blank">Test Database (110.7 Mb)</a>
         </p>

   .. container:: coltwo

      .. raw:: html

         <p>
            <a class="reference internal" href="https://doi.org/10.26180/19653840" target="_blank">expam RefSeq (122.35 Gb)</a>
         </p>

* Unzip the compressed file. For instance, if you downloaded :code:`Test Database` (:code:`test.tar.gz`) into your home directory :code:`~`, run

.. code-block:: console

    $ tar -xvzf test.tar.gz

The database will now be located at :code:`~/test`, which is the directory you should pass to any :code:`-db` flag as input to **expam**.

* Download these metagenomic reads (simulated reads from a genome used to build :code:`Test Database`).

.. container:: wideflexcontainer

   .. container:: colone

      .. raw:: html

         <p>
            <a class="reference internal" href="https://drive.google.com/file/d/1hTOndUelxf1cEEW8EIRYZrdxaHKez9qz/view?usp=sharing" target="_blank">Fasta reads (432 Kb)</a>
         </p>

* Unzip these metagenomic reads into a new folder, which we will call :code:`reads`. Assuming you have downloaded and moved the above reads into your home directory :code:`~`, run

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

* Although **expam** has the capacity to convert phylogenetic classifications into taxonomy, the default behaviour is only to provide phylogenetic results. The next section details how to obtain a taxonomic summary.

* :code:`expam classify` will create the following folder structure, for a classification run using the flag :code:`expam classify ... --out my_run`:

 ============================== =============================== ========================================================================== =================================================================================================================================================================================== 
 Folder/File                    Name                            Description                                                                Columns                                                                                                                                                                            
 ============================== =============================== ========================================================================== =================================================================================================================================================================================== 
 **my_run/phy/**                Phylogenetic output             Phylogenetic output parent folder                                          N/A                                                                                                                                                                                
 my_run/phy/classified.csv      Classifications summary file    Raw classification counts at phylogenetic nodes for each sample.           One column for each input sample.                                                                                                                                                  
 my_run/phy/split.csv           Splits summary file             Split classification counts at each phylogenetic node                      One column for each input sample.                                                                                                                                                  
 my_run/phy/[SAMPLE].csv        Sample summary file             Overview of phylogenetic classification results for a particular sample    Cumulative Classified Percentage, Cumulative Classified Count, Raw Classified Count, Cumulative Split Percentage, Cumulative Split Count, Raw Split Count                          
 **my_run/phy/raw**             Phylogenetic read-wise output   Phylogenetic read-wise output parent folder                                N/A                                                                                                                                                                                
 my_run/phy/raw/[SAMPLE].csv    Phylogenetic read-wise output   Phylogenetic read-wise output for a particular sample                      Classification Code, Read ID, Assigned Phylogenetic Node, Total Read Length, K-mer Distribution.                                                                                   
 ============================== =============================== ========================================================================== =================================================================================================================================================================================== 

* When :code:`expam classify` is run, the :code:`phy` folders are populated with results.

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
    Loading reads from reads/GCF_000005845.2_ASM584v2_genomic.fna.gz_2.fa, reads/GCF_000005845.2_ASM584v2_genomic.fna.gz_1.fa...
    The ETE3 package is not installed. Either install this module, or use the
        --itol flag if you wish to use our native iTOL integration.
    Skipping plotting of my_run/phylotree_classified.pdf...
    Skipping plotting of my_run/phy/split.csv - no samples with counts in this matrix.

.. note::

    Note that **expam** tried to plot the results on a phylotree, but since we did not have the ete3 module installed,
    it simply skipped plotting the results. This is the expected behaviour to let you know **expam** was not able
    to produce a graphical picture for your results.

    It is also possible to visualise the results with iTOL instead of the default ete3 module. To use this 
    capacity, simply supply the :code:`--itol` flag with the :code:`expam classify` command and *expam* will
    output files that serve as input to the iTOL web portal. See :ref:`itol integration`.


* The phylogenetic classifications will be located at :code:`~/my_run/phy`, and will contain four files:
    * :code:`~/my_run/phy/GCF_000005845.2_ASM584v2_genomic.gz_1.csv` - sample summary file,

    =================================== =================================== ============================== ======================= ============================== ========================= ================== 
    Node                                Cumulative Classified Percentage    Cumulative Classified Count    Raw Classified Count    Cumulative Split Percentage    Cumulative Split Count    Raw Split Count   
    =================================== =================================== ============================== ======================= ============================== ========================= ================== 
    unclassified                        0.0%                                0                              0                       0.0%                           0                         0                 
    p1                                  100.0%                              1000                           3                       0.0%                           0                         0                 
    p2                                  99.7%                               997                            232                     0.0%                           0                         0                 
    p5                                  76.5%                               765                            0                       0.0%                           0                         0                 
    GCF_000005845.2_ASM584v2_genomic    76.5%                               765                            765                     0.0%                           0                         0                 
    =================================== =================================== ============================== ======================= ============================== ========================= ================== 


    * :code:`~/my_run/phy/classified.csv` - classified summary file,

    =================================== ======================================== 
    Node                                GCF_000005845.2_ASM584v2_genomic.gz_1   
    =================================== ======================================== 
    unclassified                        0                                       
    p1                                  3                                       
    p2                                  232                                     
    GCF_000005845.2_ASM584v2_genomic    765                                     
    =================================== ======================================== 

    * :code:`~/my_run/phy/split.csv` - split summary file,

    =================================== ======================================== 
    Node                                GCF_000005845.2_ASM584v2_genomic.gz_1   
    =================================== ======================================== 
    p1                                  0                                       
    p2                                  0                                     
    GCF_000005845.2_ASM584v2_genomic    0                                     
    =================================== ======================================== 

    * :code:`~/my_run/phy/raw` - raw read-wise classifications. There will be a single raw read-wise output file, :code:`~/my_run/phy/raw/GCF_000005845.2_ASM584v2_genomic.gz_1.csv`.

    .. code-block::

        C       R6024166953890296505    p2      302     p2:240
        C       R5637238631728068726    p10     302     p10:33 p2:174 p10:33
        C       R4776396741200842676    p10     302     p10:240
        C       R5978962780799406918    p10     302     p10:240
        C       R5054328572910484091    p2      302     p2:240
        C       R3752077777745312170    p2      302     p2:240

The sample summary file is a comma-separated document where the first element of each row is a phylogenetic node/clade, and the corresponding values contain details of the raw and cumulative classifications and splits at this particular node.

The classified summary file is a comma-separated matrix where each row is a phylogenetic clade, each column is an input sample, and the cell value is the raw counts at this clade. The split summary file is an analogous file that contains the raw split count at any given clade. These two files are formatted such that they will always have the same column and row indices, and in the same order.

The raw read-wise output is a sub-directory containing one output file for each input sample, outlining read-wise output in kraken format.

A more comprehensive overview is given :doc:`this tutorial <tutorials/classify>`.


Convert to taxonomy
-------------------

* First run :code:`expam download_taxonomy` to download the taxonomy for all sequences in the database. This will require an internet connection.

.. code-block:: console

    $ expam download_taxonomy -db ~/test
    Posting 6 UIDs to NCBI Entrez nuccore.
    Received 6 response(s) for ESummary TaxID request!
    Posting 6 UIDs to NCBI Entrez taxonomy.
    Received 6 response(s) for EFetch Taxon request!
    Taxonomic lineages written to ~/test/phylogeny/taxid_lineage.csv!
    Taxonomic ranks written to ~/test/phylogeny/taxa_rank.csv!

* The taxonomic results folder has a similar format to the corresponding phylogenetic output. Again using the example with output named :code:`my_run`, **expam** will create the following folder structure.

============================== =============================== ========================================================================== =================================================================================================================================================================================== 
Folder/File                    Name                            Description                                                                Columns                                                                                                                                                                            
============================== =============================== ========================================================================== =================================================================================================================================================================================== 
**my_run/tax/**                Taxonomic output                Taxonomic output parent folder                                             N/A                                                                                                                                                                                
my_run/tax/[SAMPLE].csv        Taxonomic sample summary        Overview of taxonomic classification results for a particular sample       Cumulative Classified Percentage, Cumulative Classified Count, Raw Classified Count, Cumulative Split Percentage, Cumulative Split Count, Raw Split Count, Rank, Scientific Name   
**my_run/tax/raw**             Taxonomy read-wise output       Taxonomic read-wise output parent folder                                   N/A                                                                                                                                                                                
my_run/tax/raw/[SAMPLE].csv    Taxonomic read-wise output      Taxonomic read-wise output for a particular sample                         Classification Code, Read ID, Assigned Taxonomic ID, Read Length                                                                                                                   
============================== =============================== ========================================================================== =================================================================================================================================================================================== 

* We saved our previous classification results to :code:`~/my_run`. This is the directory we pass to :code:`expam to_taxonomy` to convert phylogenetic classifications to taxonomy.

.. code-block:: console

    $ expam to_taxonomy -db ~/test --out ~/my_run
    * Initialising node pool...
    * Checking for polytomies...
        Polytomy (degree=3) detected! Resolving...
    * Finalising index...


* There will now be taxonomic output files located in :code:`~/my_run/tax/`, analogous to each of the files present in the phylogenetic output, with the exception of :code:`classified.tsv` and :code:`split.tsv` - only the sample summaries and raw read-wise output are converted.

* :code:`~/my_run/tax/GCF_000005845.2_ASM584v2_genomic.gz_1.csv` - taxonomic summary file

=============== =================================== ============================== ======================= ============================== ========================= ================== =============== ================================================================= 
Node            Cumulative Classified Percentage    Cumulative Classified Count    Raw Classified Count    Cumulative Split Percentage    Cumulative Split Count    Raw Split Count    Rank            Scientific Name                                                  
=============== =================================== ============================== ======================= ============================== ========================= ================== =============== ================================================================= 
unclassified    0.0%                                0                              0                       0.0%                           0                         0                  0               0                                                                
1               100.0%                              1000                           0                       0.0%                           0                         0                  root                                                                             
131567          100.0%                              1000                           0                       0.0%                           0                         0                  top             cellular organisms                                               
2               100.0%                              1000                           235                     0.0%                           0                         0                  superkingdom    cellular organisms Bacteria                                      
1224            76.5%                               765                            0                       0.0%                           0                         0                  phylum          cellular organisms Bacteria Proteobacteria                       
1236            76.5%                               765                            0                       0.0%                           0                         0                  class           cellular organisms Bacteria Proteobacteria Gammaproteobacteria   
=============== =================================== ============================== ======================= ============================== ========================= ================== =============== ================================================================= 

* :code:`~/my_run/tax/raw/GCF_000005845.2_ASM584v2_genomic.gz_1.csv` - taxonomic read-wise output. The second column is the read header, the third is the assigned taxid, and the fourth is the length of the read. Observe length of 302 for paired-end 150bp reads (reads are concatenated with 'N's).

.. code-block::

    C       R6024166953890296505    2       302
    C       R5637238631728068726    511145  302
    C       R4776396741200842676    511145  302
    C       R5978962780799406918    511145  302
    C       R5054328572910484091    2       302
    C       R3752077777745312170    2       302
    C       R3409307205017145485    511145  302
    C       R6600509248827337399    2       302
    C       R9030130456493509712    511145  302

The complete comprehensive overview is given :doc:`this tutorial <tutorials/classify>`.

