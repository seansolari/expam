import io
import os
import re
from math import sqrt
from multiprocessing.shared_memory import SharedMemory
from os import ftruncate

import numpy as np
import tables as tb
from tables import UInt64Atom, UInt16Atom
from expam.sequences import format_name, get_opener

CHUNK_SIZE = int(1e9)  # 1Gb sequence file chunking.


class MaskStore:
    def __init__(self):
        self.store = {}

    def __setitem__(self, key, value):
        # Place value in store.
        self.store[key] = value

    def __getitem__(self, key):
        return self.store[key]

    def __delitem__(self, key):
        del self.store[key]


class SequenceStore:
    def __init__(self, k, out_dir):
        self._store = {}
        self.stream = io.BytesIO()
        self.k = k
        self.out_dir = out_dir

        self.accession_ids = {}

    def add_sequence(self, file_path):
        # Check it is a valid file.
        if not os.path.isfile(file_path):
            raise ValueError("Directory %s is not a file!" % file_path)

        print("Extracting sequences from %s..." % file_path)

        # Determine opening mode based on compression of sequence file.
        file_name = format_name(file_path)
        file_name, (mode, opener) = get_opener(file_name)

        with opener(file_path, mode=mode) as f:
            # Extract the accession id.
            metadata = f.readline()
            assert b'>' in metadata

            accession_id, tax_id = self.extract_accession_id(metadata)

            if accession_id is None and tax_id is None:
                raise ValueError("Couldn't identify sequence %s!" % file_name)

            for line in f:
                line = line.strip()
                # Skip information line at start of FASTA format.
                if b'>' in line:
                    self.stream.write(b'N')
                else:
                    self.stream.write(line)

        # Convert to numeric sequence and put in store.
        sequence_view = self.stream.getvalue()
        self._store[file_name] = sequence_view

        # Save accession ID for this file.
        self.accession_ids[file_name] = [accession_id, tax_id]

        # Clear stream.
        self.stream.seek(0)
        self.stream.truncate()

        return file_name, len(self._store[file_name])

    def __getitem__(self, key):
        # If it is a branch, this is a valid request,
        # but no information from the Sequence Store is required.
        if key.isdigit():
            return None
        else:
            return self._store.__getitem__(key)

    def __setitem__(self, key, value):
        self._store.__setitem__(key, value)

    def __delitem__(self, key):
        self._store.__delitem__(key)

    def keys(self):
        return self._store.keys()

    @staticmethod
    def extract_accession_id(string):
        def try_convert(byte_string):
            try:
                return byte_string.decode("utf8")

            except AttributeError:
                return byte_string

        accn_id = re.findall(b'(?:^\s*\>)(.*?)(?:\s+|\|)', string)
        tax_id = re.findall(b'(?:\|kraken:taxid\||\|taxid\|)(\d+)(?:\s+|\|)', string)

        accn_id = None if not accn_id else try_convert(accn_id[0])
        tax_id = None if not tax_id else try_convert(tax_id[0])

        return accn_id, tax_id

    def save_accession_ids(self, suffix=""):
        id_data = "\n".join([
            "%s,%s,%s" % (name, keys[0], keys[1])
            for name, keys in self.accession_ids.items()
        ])

        with open(os.path.join(self.out_dir, "phylogeny", "accession_ids_%s.csv" % suffix), "w") as f:
            f.write(id_data)


class ResizableSharedMemory(SharedMemory):
    def __init__(self, name=None, create=False, size=0):
        super(ResizableSharedMemory, self).__init__(name, create, size)

    @property
    def size(self):
        """Size in bytes."""
        return self._size

    @size.setter
    def size(self, val):
        self._size = val

    def resize(self, newsize):
        fd = self._fd

        try:
            ftruncate(fd, newsize)
            return 1, newsize
        except OSError as e:
            return 0, e

    def remap(self, newsize):
        if newsize != self.size:
            del self._buf  # Clear buffer exposure.
            self._mmap.resize(newsize)  # Remap memory.
            self._buf = memoryview(self._mmap)  # Re-expose buffer.
            # Update size.
            self.size = newsize


class TablesDb:
    def __init__(self, db_url, shape=None, keys_type=None, values_type=None, create=False, edit=False):
        self.url = db_url

        self.keys_shape = shape
        self.values_shape = None if shape is None else tuple((shape[0],))

        self.keys_type = keys_type
        self.values_type = values_type

        self.edit = True if (create or edit) else False

        self.keys_offset = 0
        self.values_offset = 0

        filters = tb.Filters(complevel=5, complib="zlib")

        if create:
            if any([shape is None, keys_type is None, values_type is None]):
                raise ValueError("Parameters not provided for creation!")

            self.fh = tb.open_file(db_url, mode="w")

            self.keys = self.fh.create_carray(self.fh.root, "keys", UInt64Atom(),
                                              self.keys_shape, filters=filters)
            self.values = self.fh.create_carray(self.fh.root, "values", UInt16Atom(),
                                                self.values_shape, filters=filters)

        else:
            if edit:
                self.fh = tb.open_file(db_url, mode="r+")

                self.keys = self.fh.root.keys
                self.values = self.fh.root.values

            else:
                self.fh = tb.open_file(db_url, mode="r")

    def _write(self, src, dest, offset):
        if not self.edit:
            raise IOError("Opened in read mode!")

        src_length = src.shape[0]
        dest[offset:offset + src_length] = src[:]

        return src_length

    def write_keys(self, arr, offset=0):
        return self._write(src=arr, dest=self.keys, offset=offset)

    def write_values(self, arr, offset=0):
        return self._write(src=arr, dest=self.values, offset=offset)

    def read_keys(self, arr):
        if self.edit:
            raise IOError("Opened in create mode!")

        self.fh.root.keys.read(out=arr)

    def read_values(self, arr):
        if self.edit:
            raise IOError("Opened in create mode!")

        self.fh.root.values.read(out=arr)

    def close(self):
        self.fh.close()


class LTMatrix:
    def __init__(self, n):
        self._width = n-2
        self._len = int((n-2) * (n-1) / 2)
        self._arr = np.zeros(self._len, dtype=np.uint16)

    def __getitem__(self, index):
        if len(index) != 2:
            raise IndexError("Can only take indices of length 2!")

        i, j = index
        if i == 0 or j == 0:
            return 0

        if i == j:
            return i

        return self._arr[self._calc_index(i - 1, j - 1)]

    def __setitem__(self, index, value):
        if len(index) != 2:
            raise IndexError("Can only take indices of length 2!")

        i, j = index
        if (i == j) or any((i == 0, j == 0)):
            raise IndexError("Can only edit lower triangular numbers!")

        self._arr[self._calc_index(i - 1, j - 1)] = value

    def __str__(self):
        string = ""

        i = 0
        for a in range(1, self._width + 1):
            string += "\t".join(str(j) for j in self._arr[i:i+a]) + "\n"
            i += a

        return string.strip()

    @staticmethod
    def _calc_index(i, j):
        if i < j:
            raise IndexError("Not lower triangular index!")

        return int((0.5 * (i-1) * i) + j)

    def save(self, out_dir):
        np.save(out_dir, self._arr)

    @classmethod
    def load(cls, in_dir):
        arr = np.load(in_dir)
        l = len(arr)
        n = int(0.5 * (3 + sqrt(1 + (8*l))))

        obj = cls(n)
        obj._arr[:] = arr[:]

        return obj
