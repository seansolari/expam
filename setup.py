"""

███████╗██╗░░██╗██████╗░░█████╗░███╗░░░███╗
██╔════╝╚██╗██╔╝██╔══██╗██╔══██╗████╗░████║
█████╗░░░╚███╔╝░██████╔╝███████║██╔████╔██║
██╔══╝░░░██╔██╗░██╔═══╝░██╔══██║██║╚██╔╝██║
███████╗██╔╝╚██╗██║░░░░░██║░░██║██║░╚═╝░██║
╚══════╝╚═╝░░╚═╝╚═╝░░░░░╚═╝░░╚═╝╚═╝░░░░░╚═╝


Metagenomic analysis using a reference phylogeny.


Uploading to PyPi
=================

Run ./buildandtwine.sh Shell script to build source code and upload.


Uploading to bioconda
=====================

Create skeleton recipe from PyPi repo.
    conda skeleton pypi --extra-specs numpy --extra-specs Cython expam


"""

import os
from setuptools import setup, find_packages
from setuptools.extension import Extension
import sys

from Cython.Build import cythonize
import numpy as np

SOURCE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SOURCE)

from src.expam import __version__


# Get project description.
with open(os.path.join(SOURCE, "README.md"), mode="r", encoding="utf-8") as f:
    long_description = f.read()

# Extension instances for Cython scripts.
extensions = [
    *cythonize(
        Extension(
            "kmers._build.kmers",
            sources=["src/expam/ext/kmers/extract.pyx", "src/expam/ext/kmers/kmers.c", "src/expam/ext/kmers/jellyfish.c"],
            include_dirs=[np.get_include()],
            extra_compile_args=["-std=c99"],
            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
        ),
        language_level="3",
        build_dir="src/expam/ext/kmers/_build"
    ),
    *cythonize(
        Extension(
            "map._build.map",
            sources=["src/expam/ext/map/map.pyx"],
            include_dirs=[np.get_include()],
            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
        ),
        language_level="3",
        build_dir="src/expam/ext/map/_build"
    ),
    *cythonize(
        Extension(
            "sets._build.sets",
            sources=["src/expam/ext/sets/sets.pyx", "src/expam/ext/sets/mfil.c"],
            include_dirs=[np.get_include()],
            extra_compile_args=["-std=c99"],
            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
        ),
        language_level="3",
        build_dir="src/expam/ext/sets/_build"
    )
]

setup(
    #
    # Metadata.
    #
    name="expam",
    version="%d.%d.%d" % __version__,
    description="Metagenomic profiling using a reference phylogeny",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/seansolari/expam",
    author="Sean Solari",
    author_email="sean.solari@monash.edu",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    #
    # Source files.
    #
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    #
    # External dependencies.
    #
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.22.0",
        "matplotlib",
        "pandas",
        "psutil",
        "requests",
        "tables"
    ],
    extras_require={
        "treeplot": ["ete3", "PyQt5==5.12", "sourmash"]
    },
    #
    # Cython modules.
    #
    ext_package="expam.ext",
    ext_modules=extensions,
    #
    # Make main callable from console.
    #
    entry_points={
        "console_scripts": [
            "expam=expam.main:main",
            "expam_limit=expam.sandbox.main:main"
        ],
    },
    #
    # Source links.
    #
    project_urls={
        "Bug Reports": "https://github.com/seansolari/expam/issues",
    },
)
