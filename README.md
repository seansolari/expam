![expam logo](docs/source/expam-logo.png)

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
user@computer:~$ git clone git@github.com:seansolari/expam.git
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

See the Quick Start Tutorial for a guide to expam's basic usage and download links for pre-built databases.

[Quick Start Tutorial](https://expam.readthedocs.io/en/latest/quickstart.html)


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

### OOM Killer

If you run into the unlikely circumstance where the OOM killer has been invoked and the program experiences an ungraceful exit, the operating system may not have cleaned all shared memory resources expam used, leading to potentially problematic memory leaks.

To prevent this occurring, make prudent use of the `expam_limit` functionality (see documentation), and don't use an extremely high number of processes (particularly for large databases). Within the range of 10-30 processes will likely be suitable for high-memory machines.

If you suspect that OOM killer has been invoked, this can be confirmed using the following command:

```bash
dmesg -T | egrep -i 'killed process'
```

In the event OOM killer has been called, it is prudent to check
how much shared memory is currently being used by the system.

```bash
df -h /dev/shm
```

If the amount of shared memory used is higher than you would expect, you can first check if there are any residual resources that need to be cleaned up.

```bash
ls -lah /dev/shm
```

If there are files starting with 'psm' and owned by you, these may be residual files that need to be cleaned up. Contact your systems administrator to remove these files.

It may also be the case that OOM killer has killed some child process, leaving the parent process sleeping (and therefore holding onto resources). You will need your system administrator's assistance to clean this up. 

To check for sleeping (expam) processes, run 

```bash
sudo lsof /dev/shm | grep "expam"
```

These sleeping processes should then be killed by running

```bash
kill -9 <PID>
```

Confirm that the leaked memory has been freed by running `df -h /dev/shm`.


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
