.. expam documentation master file, created by
   sphinx-quickstart on Thu Feb 24 09:37:07 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: expam-logo.png
    :width: 500
    :align: center
    :alt: logo

Welcome to the **expam** documentation!
=======================================

.. rst-class:: page-break

**expam** is a Python package for phylogenetic analysis of metagenomic data.

.. include:: .special.rst

.. container:: flexcontainer

   .. container:: colone

      :ref:`Installation Instructions`

   .. container:: coltwo

      :doc:`Quickstart <quickstart>`

   .. container:: colthree

      :doc:`Tutorials <tutorials/index>`

   .. container:: colfour

      .. raw:: html

         <p>
            <a class="reference internal" href="https://github.com/seansolari/expam" target="_blank">GitHub</a>
         </p>

   .. container:: colfive

      .. raw:: html

         <p>
            <a class="reference internal" href="https://figshare.com/s/3475c3a9aa926a40c722" target="_blank">Database</a>
         </p>

   .. container:: colsix

      .. raw:: html

         <p>
            <a class="reference internal" href="https://github.com/seansolari/expam/issues" target="_blank">Report Bug</a>
         </p>


.. _Installation Instructions:

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


From GitHub source
^^^^^^^^^^^^^^^^^^

To install from source, you need a local installation of Python >=3.8, as well as *numpy* and *cython*.
There are some commonly encountered problems when installing on Linux, the most common of which are
outlined in the FAQ section on the `GitHub page <https://github.com/seansolari/expam>`_.

First download the source code from the GitLab repository.

.. code-block:: console

   $ git clone git@github.com:seansolari/expam.git

This can then be installed locally by executing the following command from the source code root.

.. code-block:: console

   $ python3 setup.py install


Usage
-----

expam's CLI uses the same structure for all commands and operations:

.. code-block:: console

   $ expam [command] [args...]


For a comprehensive list of commands and arguments, see :doc:`commands <commands>`. Practical usage
of these commands for building and classifying are given in the :doc:`tutorials <tutorials/index>`.

.. image:: expam-figure-v4.2.jpg
    :align: center
    :alt: expam Pipeline.


**Important** - monitoring memory usage
---------------------------------------

Be aware of the built-in tools for monitoring and restricting **expam**'s memory usage,
outlined :ref:`here <limiting resource usage>`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart
   Documentation<commands>
   dependencies
   tutorials/index
   tree

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
