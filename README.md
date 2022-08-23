![expam logo](docs/source/expam-logo.png)

## **Install**.

#### From Bioconda (Recommended)

Conda installation is recommended, and best practise is to install *expam* in a new environment. Some users may wish to use the ETE3 toolkit for plotting, while others may prefer the iTOL tool. Both commands are included in respective order.

**With ETE3**
```console
conda create -n expam -c conda-forge -c bioconda -c etetoolkit expam ete3
```

**Without ETE3**
```console
conda create -n expam -c conda-forge -c bioconda expam
```

#### From PyPI

***Mac***

You will need a local installation of HDF5. This may already be installed on your machine, but can be installed using Homebrew with the following commands.

```console
brew install pkg-config
brew install hdf5
```

If you encounter any errors, check the FAQ section on GitHub for solutions.

Then **upgrade pip** and install expam.

```console
python3 -m pip install --upgrade pip
python3 -m pip install expam
```

***Linux***

You may need to update g++ resources on your local machine. For linux, you can run the following.

```console
apt update
apt-get install build-essential
```

Then **upgrade pip** and install expam.

```console
python3 -m pip install --upgrade pip
python3 -m pip install expam
```

#### From GitLab source

To install from source, you need a local installation of `Python >=3.8`, as well as `numpy`
and `cython`. There are some commonly encountered problems when installing on Linux, the
most common of which are outlined in the FAQ section below.

You may need to update g++ resources on your local machine. For linux, you can run the following.

```console
apt update
apt-get install build-essential
```

First download the source code from the GitLab repository.
```console
git clone https://github.com/seansolari/expam.git

```

This can then be installed locally by executing the following command from the
source code root:
```console
cd expam
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
python3 setup.py install
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
```console
sudo apt-get install build-essential
```

<hr>

**fatal error: Python.h: No such file or directory**

This simply means you need to install/update the Python development files for version 3.
```console
sudo apt-get install python3-dev
```

(Reference - [SO](https://stackoverflow.com/questions/21530577/fatal-error-python-h-no-such-file-or-directory/21530768))

<hr>

**ERROR:: Could not find a local HDF5 installation** (Mac)

Ensure you have HDF5 installed using Homebrew:

```console
brew install pkg-config
brew install hdf5
```

If you see
```console
You may need to explicitly state where your local HDF5 headers and library can be found by setting the ''HDF5_DIR'' environment variable or by using the ''--hdf5'' command-line option.
```
you will need to explicitly set the `HDF5_DIR` environment variable. To see where HDF5 has been installed, run
```console
brew info hdf5
```
You should see something like `/usr/local/Cellar/hdf5/VERSION...` or `/opt/local/Cellar/hdf5/VERSION...` (ie. ignore everything after the complete version, which will have numbers separated by dots). Then set this environment variable with
```console
HDF5_DIR=/opt/local/Cellar/hdf5/VERSION,
```
replacing this path with your output from `brew info`.

Now retry the installation having set this environment variable.

<hr>

**ete3 importing errors**

For instance, `ImportError: cannot import name 'NodeStyle'`.

The *ete3* module depends on Qt to draw things, and there are two stages to getting this to work: first, Qt needs to be installed, and then you need to let Python know that Qt is installed. Follow the following instructions depending on your OS.

***Mac***

Install qt5 using brew.

```bash
brew install qt5
brew list --versions qt5
```

This should show you the precise version that brew installed. We now tell Python which version of Qt5 to link up with. Say we have `qt@5 5.15.3` from the above command, then we would run

```bash
python3 -m pip install pyqt5==5.15
```

Had the output been `qt@5 5.12.0`, we would run 

```bash
python3 -m pip install pyqt5==5.12
```

ie. the first two parts of the version from brew. This should remedy the problem.


***Linux***

First update the local installation of Qt.
```bash
sudo apt-get install qt5-default
```

Now double-check which version of Qt has been installed.
```bash
dpkg -l | grep "pyqt5"
```

Take the first two parts of the version output from this, and pass it to this following install with Pip. For instance, if we have `qt5 5.12.0`, take the `5.12` component. Install the corresponding Python interface to Qt.
```bash
python3 -m pip install pyqt5==5.12
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
expam --version
expam --help
...
```

<hr style="border:1px solid #ADD8E6"> </hr>

## Bug Reports
Please raise any bug reports at https://github.com/seansolari/expam/issues
accompanied by any error messages, a rough description of the database/setup and
parameters used to create the database.
