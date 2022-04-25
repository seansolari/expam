from multiprocessing import shared_memory

import numpy as np
import tables as tb

from expam.database import DataTypeConfig, FileLocationConfig, expam_dtypes
from expam.database.config import load_database_config


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

            self.fh: tb.File = tb.open_file(db_url, mode="w")

            self.keys = self.fh.create_carray(self.fh.root, "keys", tb.UInt64Atom(),
                                              self.keys_shape, filters=filters)
            self.values = self.fh.create_carray(self.fh.root, "values", tb.UInt16Atom(),
                                                self.values_shape, filters=filters)

        else:
            if edit:
                self.fh: tb.File = tb.open_file(db_url, mode="r+")

                self.keys = self.fh.root.keys
                self.values = self.fh.root.values

            else:
                self.fh: tb.File = tb.open_file(db_url, mode="r")

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


class SharedDB:
    def __init__(self, out_dir, keys_shape, values_shape, dtypes: DataTypeConfig = expam_dtypes, config: FileLocationConfig = None):
        print("Preparing shared memory allocations...")

        # Prepare shm segment for keys.
        keys_shm = shared_memory.SharedMemory(
            size=np.product(keys_shape) * dtypes.keys_dtype_size,
            create=True
        )
        keys_interface = np.ndarray(
            shape=keys_shape,
            dtype=dtypes.keys_dtype,
            buffer=keys_shm.buf
        )

        # Prepare shm segment for values.
        values_shm = shared_memory.SharedMemory(
            size=np.product(values_shape) * dtypes.map_dtype_size,
            create=True
        )
        values_interface = np.ndarray(
            shape=values_shape,
            dtype=dtypes.map_dtype,
            buffer=values_shm.buf
        )

        if config is None:
            config = load_database_config(out_dir)

        db = TablesDb(config.database_file)

        #
        print("Loading database keys...")
        db.read_keys(arr=keys_interface)

        #
        print("Loading database values...")
        db.read_values(arr=values_interface)

        db.close()

        # Remember the names for attaching in future.
        self.keys_shm_name = keys_shm.name
        self.values_shm_name = values_shm.name

        # Remember array data.
        self.keys_meta = {
            "shape": keys_shape,
            "dtype": dtypes.keys_dtype
        }
        self.values_meta = {
            "shape": values_shape,
            "dtype": dtypes.map_dtype
        }

        # Close my access.
        keys_shm.close()
        values_shm.close()

