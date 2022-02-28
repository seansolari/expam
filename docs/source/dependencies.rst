Additional Dependencies (optional)
==================================

Setting up a conda environment
------------------------------
Conda is an incredibly convenient tool to keep all packages and tools used by expam
in the same place, isolating it from all other tools on your computer.

See `here <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ for
instructions on how to install Conda. Conda searches online repositories for the software
you ask it to install. There are certain repositories you should let Conda know about, which
contain a lot of useful bioinformatics software.

.. code-block:: console

    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge


We now create the conda environment, and call it **expam**.

.. code-block:: console

    $ conda create --name expam python=3.8


You enter the environment using :code:`conda activate`, at which point your console will indicate
you are in an environment by adding (expam) to the start of each line. All applications
installed using Conda from within this environment will only be accessible upon entering
the environment. You can leave the environment using :code:`conda deactivate`.

.. code-block:: console

    $ conda activate expam
    (expam) $ conda list

    # packages in environment at /opt/anaconda3/envs/test:
    #
    # Name                    Version                   Build  Channel
    ca-certificates           2021.10.8            h033912b_0    conda-forge
    libcxx                    12.0.1               habf9029_0    conda-forge
    libffi                    3.4.2                h0d85af4_5    conda-forge
    libzlib                   1.2.11            h9173be1_1013    conda-forge
    ncurses                   6.2                  h2e338ed_4    conda-forge
    openssl                   3.0.0                h0d85af4_2    conda-forge
    pip                       21.3.1             pyhd8ed1ab_0    conda-forge
    python                    3.8.12          h43ca1e7_2_cpython    conda-forge
    python_abi                3.8                      2_cp38    conda-forge
    readline                  8.1                  h05e3726_0    conda-forge
    setuptools                59.1.1           py38h50d1736_0    conda-forge
    sqlite                    3.36.0               h23a322b_2    conda-forge
    tk                        8.6.11               h5dbffcc_1    conda-forge
    wheel                     0.37.0             pyhd8ed1ab_1    conda-forge
    xz                        5.2.5                haf1e3a3_1    conda-forge
    zlib                      1.2.11            h9173be1_1013    conda-forge
    (expam) $ conda deactivate
    $  


Additional Software
-------------------
If you wish to build your own tree for a custom collection of genomes using the **expam**
interface, you will need some tools. See :doc:`Tutorial Three <tutorials/tree>` for a
tutorial on how to build trees using the **expam** interface, this should make you aware
of which software you will need. Brief instructions for installing all possible
software are shown below.

Note: whenever the tutorials mention a 'local installation' of some tool, this simply
means that the tool should be available from a call in the console. Installing any required
tools within the environment using the following commands will satisfy this.

Installing Mashtree
^^^^^^^^^^^^^^^^^^^
`Mashtree GitHub <https://github.com/lskatz/mashtree>`_

Mashtree requires a local installation of mash. Run the following commands in succession.

.. code-block:: console

    $ conda install mash
    $ conda install mashtree


Installing Mash
^^^^^^^^^^^^^^^
`Mash GitHub <https://github.com/marbl/Mash>`_

.. code-block:: console

    $ conda install mash


Installing RapidNJ
^^^^^^^^^^^^^^^^^^
`RapidNJ homepage <https://birc.au.dk/software/rapidnj/>`_

.. code-block:: console

    $ conda install rapidnj


Installing QuickTree
^^^^^^^^^^^^^^^^^^^^
`QuickTree GitHub <https://github.com/khowe/quicktree>`_

.. code-block:: console

    $ conda install quicktree


