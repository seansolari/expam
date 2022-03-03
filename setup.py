import os
from setuptools import setup, find_packages
from setuptools.extension import Extension

from Cython.Build import cythonize
import numpy as np

EXPAM_VERSION = (0, 0, 8)

SOURCE = os.path.dirname(os.path.abspath(__file__))

# Get project description.
with open(os.path.join(SOURCE, "README.md"), mode="r", encoding="utf-8") as f:
    long_description = f.read()

# Extension instances for Cython scripts.
extensions = [
    Extension(
        "map",
        sources=["src/expam/c/map.pyx"],
        include_dirs=[np.get_include()],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
    ),
    Extension(
        "extract",
        sources=["src/expam/c/extract.pyx", "src/expam/c/kmers.c", "src/expam/c/jellyfish.c"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-std=c99"],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
    ),
    Extension(
        "sets",
        sources=["src/expam/c/sets.pyx", "src/expam/c/mfil.c"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-std=c99"],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
    ),
]

setup(
    #
    # Metadata.
    #
    name="expam",
    version="%d.%d.%d" % EXPAM_VERSION,
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
        # "ete3",
        # "PyQt5==5.12",
        "numpy",
        "matplotlib",
        "pandas",
        "psutil",
        "requests",
        # "sourmash",
        "tables"
    ],
    #
    # Cython modules.
    #
    ext_package="expam.c",
    ext_modules=cythonize(extensions, language_level="3"),
    #
    # Make main callable from console.
    #
    entry_points={
        "console_scripts": [
            "expam=expam.cli:main",
            "expamlimit=expam.sandbox:main"
        ],
    },
    #
    # Source links.
    #
    project_urls={
        "Bug Reports": "https://github.com/seansolari/expam/issues",
    },
)
