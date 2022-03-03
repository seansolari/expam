Tree building with **expam**
============================

**Figure 1:** Overview of tree building with expam.

.. mermaid:: treebuildingpipeline.mmd


Part One: Building a tree
-------------------------

.. note:: 

    Instructions for installing all software mentioned in this section can be found
    in :doc:`Dependencies tutorial <../dependencies>`.

* There are three stages to building a distance tree:
  1. sketching the sequences (a compressed form that still contains a lot of the important details),
  2. using these sketches to estimate sequence distance,
  3. apply Neighbour Joining (NJ) to the distance matrix to compute a distance tree.

* This process can be done manually, or automated through the :code:`mashtree` application.

Basic tree building
^^^^^^^^^^^^^^^^^^^

* We'll build a tree from the small collection of reference genomes supplied with the **expam** source code.
* First create a new database and add these sequences.
  * I will refer to the location of these sequences as :code:`../expam/test/data/sequence/`, but this may differ for you.

.. code-block:: console

    $ expam create -db tree
    $ expam add -db tree -d ../expam/test/data/sequence/

* We'll set parameters that have been used before to build the database, but setting them here will be relevant for tree building.
  
.. code-block:: console

    $ expam set -db tree -k 31 -n 4 -s 100

* This means we will be creating sequence sketches with a k-mers size :code:`k=31` and using a sketch size :code:`s=100`. 
* We will also let the tree building software use up to :code:`n=4` threads at a time.

The easy way
""""""""""""

* From here, :code:`mashtree` can automate the process for us (you must have it installed).

.. code-block:: console

    $ expam mashtree -db tree
    $ expam tree -db tree

The less easy way
"""""""""""""""""

* We'll manually go through the three steps outlined above.
* First sketch the sequences.

.. code-block:: console

    $ expam sketch -db tree
    $ ls tree/phylogeny/sketch
    default.k31.s100.msh

.. note:: 

    By default, this will use :code:`mash` to create the sketches. Alternatively, 
    provide the :code:`--sourmash` flag to use :code:`sourmash` for sketching.

* Now we create the distance matrix from these sketches.

.. code-block:: console

    $ expam distance -db tree
    $ ls tree/phylogeny/distance
    default.k31.s100.tab
    $ head -n 3 tree/phylogeny/distance/default.k31.s100.tab
    6
    GCF_000006765.1_ASM676v1_genomic.fna.gz	0	1	1	1	1	1
    GCF_000005845.2_ASM584v2_genomic.fna.gz	1	0	1	1	0.0158863	1

* This shows the tab-delimited computed distance matrix.

.. warning:: 

    If you used :code:`--sourmash` to create the sketches, you must also supply
    :code:`--sourmash` when computing the distance matrix.

* Finally, we'll use a NJ tool to compute the tree from this matrix.

.. code-block:: console

    $ expam nj -db tree

.. note:: 

    By default, **expam** relies on RapidNJ to do NJ. However, it can call a local installation
    of QuickTree using :code:`--quicktree` (if you have that installed).

    .. code-block:: console

        $ expam nj -db tree --quicktree

* The tree can now be finalised and attached to the database using the :code:`tree` command.

.. code-block:: console

    $ expam tree -db tree

.. note:: 

    Running :code:`expam sketch`, :code:`expam distance` and :code:`expam nj` is therefore
    *roughly* equivalent to :code:`expam mashtree`, at least from the perspective of outcome.


Part Two: Building a tree in parts
----------------------------------

* You may wish to build a tree containing both bacterial and viral genomes, or even some human sequences for contamination detection. These genomes are very different sizes, and so using the same :code:`k` and :code:`s` parameters for each of these types of genomes may not produce very accurate trees. It may be more prudent to build trees for each of these organism types separately, and then join these subtrees afterwards.
* **expam** implements a set of routines that enable you to construct a tree in this way - i.e. in parts.
* Say we have two groups of sequences, :code:`a` and :code:`b`, that we want to construct trees for *separately*.
* We can separate these sequences in the database by adding then to separate :code:`groups`.

.. code-block:: console

    $ expam add -db tree --group a -d ~/Documents/Sequences/genomes/a/
    $ expam add -db tree --group b -d ~/Documents/Sequences/genomes/b/

.. note:: 

    Run :code:`expam print -db tree` and notice how **expam** lists multiple groups.

.. note:: 

    Even though these sequences are added to separate groups, they are all part of the reference collection - 
    **it won't affect the later database build behaviour.**

* We will use :code:`k=31, s=1000` for group :code:`a`, and :code:`k=21, s=100` for group :code:`b`.

.. code-block:: console

    $ expam set -db tree --group a  -k 31 -s 1000
    $ expam set -db tree --group b -k 21 -s 100

.. note:: 

    By specifying parameters alongside a :code:`--group` flag, **expam** recognises that these parameters 
    are specifically for tree building, not for database construction. Those would still need to be set via

    .. code-block:: console

        $ expam set -db tree -k 31 -n 4

    If database parameters have been set (those without :code:`--group` flags) but specific group flags have 
    not been set, **expam** will automatically use the database build parameters for tree building.

* We also need to tell **expam** how these trees will be joined at the end. There are two rules for this template:
  1. It is a Newick format tree, where group names appear in double braces.
  2. The template must be placed at :code:`database_name/phylogeny/tree/database_name.nwk` (the :code:`phylogeny` subdirectory in the database folder). Replace :code:`database_name` with your database name.

* The template we will use is

.. code-block:: 

    ({{a}},{{b}});

Build in parts - mashtree
^^^^^^^^^^^^^^^^^^^^^^^^^

* We can supply the :code:`mashtree` command to build these two trees separately - **expam** takes care of this behind the scenes.

.. code-block:: console

    $ expam mashtree -db tree

* Now finalise with the :code:`tree` command to apply the template and cobine these tree.

.. code-block:: console

    $ expam tree -db tree

.. note::

    If you wanted to, you could run :code:`mashtree` on these groups separetely.

    .. code-block:: console

        $ expam mashtree --group a
        $ expam mashtree --group b

Build in parts - manual
^^^^^^^^^^^^^^^^^^^^^^^

* Despite having split the sequences into groups, running the same chain of :code:`sketch`, :code:`distance` and :code:`nj` commands will result in **expam** running these commands on each group consecutively.
* Sketch the sequences

.. code-block:: console

    $ expam sketch -db tree

which is equivalent to

.. code-block:: console

    $ expam sketch -db tree --group a
    $ expam sketch -db tree --group b

* You can confirm these two groups have sketch files.

.. code-block:: console

    $ ls tree/phylogeny/sketch
    a.k31.s1000.msh
    b.k21.s100.msh

* Get pairwise distances

.. code-block:: console

    $ expam distance -db tree

which is equivalent to

.. code-block:: console

    $ expam distance -db tree --group a
    $ expam distance -db tree --group b

* Distances can be found in the :code:`tree/phylogeny/distance/` folder.

* Finally, apply NJ

.. code-block:: console

    $ expam nj -db tree

which is equivalent to

.. code-block:: console

    $ expam nj -db tree --group a
    $ expam nj -db tree --group b

* Finalise the tree using the template.

.. code-block:: console

    $ expam tree -db tree

.. note:: 

    As in Part One, the :code:`--sourmash` and :code:`--quicktree` flags can be supplied to use
    those alternative softwares.

