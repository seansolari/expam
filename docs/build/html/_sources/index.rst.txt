.. expam documentation master file, created by
   sphinx-quickstart on Thu Feb 24 09:37:07 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to expam's documentation!
=================================

**expam** is a Python package for phylogenetic analysis of metagenomic data.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Documentation<commands>
   dependencies
   tutorials/index


Installation
------------

See the :doc:`dependencies <dependencies>` tutorial for more detailed instructions for creating
a virtual environment and managing expam's code and dependencies.

Conda
^^^^^
Conda installation is recommended.

.. code-block:: console

   $ conda install expam


Python Package Index (pip)
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

   $ pip install expam


Usage
-----

expam's CLI uses the same structure for all commands and operations:

.. code-block:: console

   $ expam [command] [args...]


For a comprehensive list of commands and arguments, see :doc:`commands <commands>`. Practical usage
of these commands for building and classifying are given in the :doc:`tutorials <tutorials/index>`.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
