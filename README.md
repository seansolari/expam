![expam logo](docs/source/expamlogo.png)

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

View our online documentation!

[https://expam.readthedocs.io/en/latest/index.html](https://expam.readthedocs.io/en/latest/index.html)


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

A complete list of available commands can by found by using the `-h`/`--help`
flags.
```console
user@computer:~$ expam --help
...
```

<hr style="border:1px solid #ADD8E6"> </hr>

## Bug Reports
Please raise any bug reports at https://github.com/seansolari/expam/issues
accompanied by any error messages, a rough description of the database/setup and
parameters used to create the database.
