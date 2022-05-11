Available commands
==================

Setting the default database
----------------------------
The EXPAM_DEFAULT_DB environment variable defines expam's default database.

.. code-block:: console

    $ expam default_db -db DB_NAME
    export EXPAM_DEFAULT_DB=PATH_TO_DB

As a shortcut, 

.. code-block:: console

    $ expam default_db -db DB_NAME >> ~/.bash_profile

sets the environment variable each time you open bash console.


Create a new database
---------------------

.. code-block:: console

    $ expam create -db DB_NAME


Set build parameters
^^^^^^^^^^^^^^^^^^^^
Set the build parameters for some database.

.. code-block:: console

    $ expam set -db DB_NAME [args...]

.. option:: -k <int>, --kmer <int>

    K-mer size for building database.

    .. note::

        A k-mer length of 31 is often used, and is probably good for most bacterial genomes.


.. option:: -n <int>, --n_processes <int>

    Number of processes to distribute building work amongst.

    .. note::

        Using more processes decreases the build time but increases the amount of computer
        resources used during build. 
        
        The following short python excerpt will tell you the **max** number of processes.

        .. code-block:: python

            >>> import multiprocessing
            >>> multiprocessing.cpu_count()
            8

.. option:: -p <file path>, --phylogeny <file path>

    Path to |newick_link| detailing tree of reference sequences.

.. |newick_link| raw:: html

    <a href="https://en.wikipedia.org/wiki/Newick_format" target="_blank">Newick file</a>


Example
"""""""

.. code-block:: console

    $ expam set -db /path/to/db -k 31 -n 12 -p /home/seansolari/tree.nwk


Add/Remove sequences
^^^^^^^^^^^^^^^^^^^^
Add reference sequences to the database.

.. code-block:: console

    $ expam add -db DB_NAME [args...]
    $ expam remove -db DB_NAME [args...]


.. option:: -d <file path>, --directory <file path>

    Add sequence at **file path** to the database.

    .. note::

        File path can be a file, or a folder.

        If a folder, expam will add all sequences within this folder.

.. option:: --first_n <int>

    *(optional)*

    Add first **n** (order same as appears in `ls`) sequences from directory into the database.

.. option:: --group <str>

    *(optional)*

    Add sequences to particular sequence group.

    See :doc:`Tutorial 1 <tutorials/overview>` for details.


Examples
""""""""

.. code-block:: console

    $ expam add -db /path/to/db -d /path/to/sequence.fna
    Added 1 sequence from /path/to/sequence.fna!

    $ expam add -db /path/to/db -d /path/to/folder/
    Added 19 sequences from /path/to/folder!


Printing database parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Print current database configuration.

.. code-block:: console

    $ expam print -db DB_NAME
    <<< expam configuration file: DB_NAME >>>

    phylogeny       -->     None
    k               -->     31
    n               -->     12
    sketch          -->     None
    pile            -->     512

    ----------------
    group name: default
            k               -->     None
            sketch          -->     None
            sequences       -->     9


Build a database
----------------
When all parameters have been set, a database can be built.

.. code-block:: console

    $ expam build -db DB_NAME


Classification
--------------

Classify
^^^^^^^^
Run metagenomic reads against a succesfully built database. See :doc:`Tutorial 2 <tutorials/classify>` for more details.

.. code-block:: console

    $ expam classify -db DB_NAME [args...]

.. option:: -d <file path>, --directory <file path>

    Path containing fasta/fastq files to be classified. 

    .. note::

        Path here can either be a single file, or a folder in which case expam
        will process all sequence files in this folder.

.. option:: --paired

    To be supplied when sample files contained paired-end reads.

.. option:: -o <str>, --out <str>

    Path to save classification results and output in.

.. option:: --taxonomy

    Convert phylogenetic results to taxonomic results.

    .. note:: 

        This requires taxonomic information for all reference sequences, see the
        :ref:`download_taxonomy <download taxonomy>` command.

.. option:: --cpm <int>, --cutoff <int>

    Apply cutoff to read count. *cpm* is short for counts-per-million, and takes priority if both *cpm*
    and *cutoff* are supplied.

.. option:: --phyla

    Colour phylotree results by phyla.

.. option:: --keep_zeros

    Keep nodes in output where no reads have been assigned.

.. option:: --ignore_names

    Don't plot names of reference genomes in output phylotree.

.. option:: --colour_list <hex string> <hex string> ...

    List of colours to use when plotting groups in phylotree.

.. option:: --group <sample name> <sample name> ...

    Space-separated list of sample files to be treated as a single group in phylotree.
    Groups are explained in this :ref:`tutorial <groups explanation>`.

    .. note::

        You may also supply a hex colour code directly after the *--group* flag
        to assign some colour code to this group of samples.

        .. code-block:: console

            $ expam classify ... --group #FF0000 sample_one sample_two

.. option:: --alpha <float>

    Percentage requirement for classification subtrees (see :doc:`Tutorial 1 <tutorials/overview>`
    and :doc:`Tutorial 2 <tutorials/classify>`).

.. option:: --itol

    Rather than use :code:`ete3` for plotting the phylogenetic tree, **expam** will output files that can be
    used with iTOL for plotting. See the :ref:`classification tutorial <itol integration>` for details.

.. option:: --log-scale

    Compute a log-transform on the counts at each node in the phylogenetic tree before 
    depiction on the phylotree.

    .. note::

        For a given sample :math:`S`, with minimum and maximum counts :math:`\underline{c}` and :math:`\overline{c}` 
        respectively (:math:`\underline{c} > 0` i.e. the smallest non-zero score), the log-transform :math:`f` of some count :math:`x` is defined by

        .. math::

            f(x) = \frac{ \log\left(x / \underline{c}\right) }{ \log\left(\overline{c} / \underline{c}\right) },

        so that :math:`f(x)\in[0,1]`. Then :math:`f(x)` is treated as an opacity score for plotting purposes.


Example
"""""""

.. code-block:: console

    $ expam classify -db DB_NAME -d /path/to/paired/reads --paired --out ~/paired_reads_analysis --taxonomy

.. _download taxonomy:

Download taxonomic data
^^^^^^^^^^^^^^^^^^^^^^^^

Download taxonomic metadata associated with all sequences used to build the database.

.. code-block:: console

    $ expam download_taxonomy -db mydb

.. note::
    This command can only be run after the database has been built. This is because
    **expam** first finds NCBI accession IDs or explicit taxon IDs in the header of each
    reference sequence, and uses these to search against the NCBI Entrez database.

.. note::
    The NCBI taxonomic ID of a reference sequence can be explicitly stated in the format

    .. code-block:: text
        
        > accessionid or metadata|taxid|TAXID|

    For instance,

    .. code-block:: text

        >NZ_RJQC00000000.1|taxid|2486576|

    If both accession ID and taxon ID are supplied, taxon ID takes precedence.

Convert results to taxonomy
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Translate phylogenetic classification output to NCBI taxonomy.

.. code-block:: console

    $ expam to_taxonomy --db DB_NAME

Plotting results on phylotree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Results are automatically visualised on top of a phylogenetic tree when during the :code:`expam classify` command,
but can also be done after classification using the :code:`phylotree` command.

.. code-block::

    $ expam phylotree -db DB_NAME --out /path/to/classification/output [args...]

.. option:: -o <str>, --out <str>

    Path to retrieve classification results for plotting.

.. option:: --phyla

    Colour phylotree results by phyla.

.. option:: --keep_zeros

    Keep nodes in output where no reads have been assigned.

.. option:: --ignore_names

    Don't plot names of reference genomes in output phylotree.

.. option:: --colour_list <hex string> <hex string> ...

    List of colours to use when plotting groups in phylotree.

.. option:: --group <sample name> <sample name> ...

    Space-separated list of sample files to be treated as a single group in phylotree.
    Groups are explained above, and in this :ref:`tutorial <groups explanation>`.

.. option:: --itol

    Rather than use :code:`ete3` for plotting the phylogenetic tree, **expam** will output files that can be
    used with iTOL for plotting. See the :ref:`classification tutorial <itol integration>` for details.

.. option:: --log-scale

    Compute a log-transform on the counts at each node in the phylogenetic tree before 
    depiction on the phylotree.

.. _limiting resource usage:

Limiting resource usage
-----------------------

**expam** allows you to provide an :code:`expam_limit` context before the :code:`expam` call to limit
how much RAM is used. *Note that this doesn't change any underlying algorithms, it simply
prepares a graceful exit of the program if it exceeds the supplied limit.* See :ref:`examples<limit example>`
for an example usage.

.. option:: -m <int>, --memory <int>

    Memory limit in bytes.

.. option:: -x <float>, --x <float>

    Percentage of total available memory to limit to.

.. option:: -t <float>, --interval <float>

    Intervals in which program memory usage is written to log file.

.. option:: -o <str>, --out <str>

    Log file to write to. By default, logs are written to console.

.. _limit example:

Example
^^^^^^^
 The following will perform a database build while restricting *expam*'s total
 memory usage to half of the available machine's RAM, writing logs
 in 1 second intervals to a :code:`build.log` file.

 .. code-block:: console

     $ expam_limit -x 0.5 -t 1.0 -o build.log expam build ...

.. warning::

    It is important that the :code:`expam_limit` command comes before
    the :code:`expam` command.

.. note::

    The :code:`expam_limit` context works the same for any command. :code:`expam build`
    can be replaced with :code:`expam classify`, or any other command.

The following is an example of the (tab-separated) log file output:

.. code-block::

    2022-03-11 02:25:05,888 ... 	total	used	free	shared	buff/cache	available
    2022-03-11 02:25:05,903 ... Mem:	944Gi	1.6Gi	427Gi	0.0Ki	515Gi	938Gi
    2022-03-11 02:25:06,915 ... Mem:	944Gi	1.6Gi	427Gi	0.0Ki	515Gi	938Gi
    2022-03-11 02:25:07,928 ... Mem:	944Gi	2.2Gi	427Gi	38Mi	515Gi	937Gi
    2022-03-11 02:25:08,940 ... Mem:	944Gi	2.2Gi	426Gi	195Mi	515Gi	937Gi
    2022-03-11 02:25:09,953 ... Mem:	944Gi	2.2Gi	426Gi	353Mi	515Gi	937Gi
    2022-03-11 02:25:10,966 ... Mem:	944Gi	2.2Gi	426Gi	516Mi	516Gi	937Gi
    2022-03-11 02:25:11,980 ... Mem:	944Gi	2.2Gi	426Gi	682Mi	516Gi	936Gi
    2022-03-11 02:25:12,992 ... Mem:	944Gi	2.2Gi	426Gi	848Mi	516Gi	936Gi
    2022-03-11 02:25:14,005 ... Mem:	944Gi	2.2Gi	425Gi	1.0Gi	516Gi	936Gi

