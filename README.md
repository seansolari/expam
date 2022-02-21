<pre style="text-align:center; background-color:rgba(0, 0, 0, 0.04); border: thin dotted black; padding:25px;">
 
    __________   ___ .______      ___      .___  ___. 
    |   ____\  \ /  / |   _  \    /   \     |   \/   | 
    |  |__   \  V  /  |  |_)  |  /  ^  \    |  \  /  | 
    |   __|   >   <   |   ___/  /  /_\  \   |  |\/|  | 
    |  |____ /  .  \  |  |     /  _____  \  |  |  |  | 
    |_______/__/ \__\ | _|    /__/     \__\ |__|  |__|

</pre>

### Contents

[**Install**.](#--install--)
- [From PyPI.](#from-pypi)
- [From GitLab source.](#from-gitlab-source)

[**Tutorials**.](#--tutorials--)

[**FAQ**](#--FAQ--)

[**Commands**.](#--commands--)

[**Bug Reports**.](#--bug-reports--)

<hr style="border:1px solid #ADD8E6"> </hr>

## **Install**.

#### From Bioconda (Recommended)

#### From PyPI

```console
user@computer:~$ pip install expam
```

#### From GitLab source

To install from source, you need a local installation of `Python >=3.8`, as well as `numpy`
and `cython`. There are some commonly encountered problems when installing on Linux, the
most common of which are outlined in the FAQ section below.

First download the source code from the GitLab repository.
```console
user@computer:~$ git clone git@gitlab.erc.monash.edu.au:ssol0002/pam.git
```
This can then be installed locally by executing the following command from the
source code root:
```console
user@computer:~$ python3 setup.py install
```

<hr style="border:1px solid #ADD8E6"> </hr>

## Documentation
[Documentation](./docs/README.md)

An outline of all available commands and flags is can be found in the `/docs/` folder.

<hr style="border:1px solid #ADD8E6"> </hr>

## Tutorials

#### 0. Setting up a conda environment (optional)
[Link](./docs/tutorials/0_setup/README.md) - setting up the conda environment and 
installing expam dependencies.

#### 1. Basics tutorial
[Link](./docs/tutorials/1_basic/README.md) - building a database and classifying reads.<br>

#### 2. Classification output in-depth
[Link](./docs/tutorials/2_classify/README.md) - organising sample classifications, and displaying results. <br>

#### 3. Building a tree
[Link](./docs/tutorials/3_tree/README.md) - create a distance tree to be used for building the database.

#### 4. Visualising results
[Link](./docs/tutorials/4_graphical/README.md) - depicting classification results on a phylotree.

<hr style="border:1px solid #ADD8E6"> </hr>

## FAQ

### Problems during installation

**error: g++: Command not found**

This is simply a matter of updating the compiler.
```bash
> sudo apt-get install build-essential
```

<hr>

**fatal error: Python.h: No such file or directory**

This simply means you need to install/update the Python development files for version 3.
```bash
> sudo apt-get install python3-dev
```

(Reference - [SO](https://stackoverflow.com/questions/21530577/fatal-error-python-h-no-such-file-or-directory/21530768))

<hr>

**ete3 importing errors**

For instance, `ImportError: cannot import name 'NodeStyle'`.

The *ete3* module depends on Qt, and for Linux it may take some tweaking to get Python
to recognise the local installation of Qt. The following seems to work for a broad
collection of circumstances.

First update the local installation of Qt.
```bash
> sudo apt-get install qt5-default
```

Now double-check which version of Qt has been installed.
```bash
> dpkl -l | grep "pyqt5"
```

Install the corresponding Python interface to Qt.
```bash
> pip3 install pyqt5==5.12
```


<hr style="border:1px solid #ADD8E6"> </hr>

## Commands

|        Flag        	|   Input   	|                                                              Description                                                             	|
|:------------------:	|:---------:	|:------------------------------------------------------------------------------------------------------------------------------------:	|
|    `-h`/`--help`   	|    N/A    	|                                                          Print help command.                                                         	|
|  `-db`/`--db_name` 	|  diretory 	|                                          Path to database directory (relative or absolute).                                          	|
|     `-k`/`--k`     	|  integer  	|                                                         K-value for analysis.                                                        	|
|     `-n`/`--n`     	|  integer  	| Number of processes spawned when building database and classifying read. <br> Number of threads passed to mashtree in `mashtree` command. 	|
| `-p`/`--phylogeny` 	| directory 	|                                                    Path to phylogeny/Newick file.                                                    	|
| `-d`/`--directory` 	| directory 	|                                                          Path to directory.                                                          	|
|    `-y`/`--pile`   	|  integer  	|         Number of genomes to pile when making database. <br> Leave as default if you don't want to mess with the natural order.        	|
|      `--first`     	|  integer  	|                                   Add first [value] genomes from some directory into the database.                                   	|


A complete list of available commands can by found by using the `-h`/`--help`
flags.
```console
user@computer:~$ expam --help
...
```

<hr style="border:1px solid #ADD8E6"> </hr>

## Bug Reports
Please raise any bug reports at https://gitlab.erc.monash.edu.au/ssol0002/pam/-/issues
accompanied by any error messages, a rough description of the database/setup and
parameters used to create the database.
