import math
import multiprocessing

import numpy as np
from expam.database import TIMEOUT, DataTypeConfig, expam_dtypes
from expam.database.db import TablesDb
from expam.ext.map import binarysearch_put
from expam.ext.sets import disjoint as make_disjoint, filler
from expam.database.shm import ResizableSharedMemory
from expam.logger import log_timer
from expam.process.jobworker import JobWorker


class UnionWorker(JobWorker):
    def __init__(
        self, main_con, lca_matrix, pile_size, logging_dir,
        lock_value, status_index, statuses, child_statuses,
        timeout=TIMEOUT, phases=None, dtype_config: DataTypeConfig = expam_dtypes
    ):
        super(UnionWorker, self).__init__(
            logging_dir, status_index, statuses, child_statuses,
            timeout=timeout, phases=phases
        )

        # Connection to main process.
        self.job_name = "union"
        self.con = main_con

        # For mapping LCAs in mapping stage.
        self.lca_matrix = lca_matrix

        # Piling parameters.
        self.in_allocations = {}
        self.in_arrays = {}
        self.pile = []

        if pile_size == "inf":
            self.pile_threshold = math.inf

        else:
            pile_size = self.is_int(pile_size)

            if pile_size is not False:
                self.pile_threshold = pile_size

            else:
                self.pile_threshold = 255  # Default.

        # Shared memory allocations.
        self.dtypes: DataTypeConfig = dtype_config
        self.shm_allocations = {
            "keys": {
                "array": None,
                "allocation": None,
                "shm_name": None,
                "size": None,
                "data_type": dtype_config.keys_dtype,
                "data_type_size": dtype_config.keys_dtype_size,
                "active": False,
                "write_offset": None,
            },
            "values": {
                "array": None,
                "allocation": None,
                "shm_name": None,
                "size": None,
                "data_type": dtype_config.map_dtype,
                "data_type_size": dtype_config.map_dtype_size,
                "active": False,
                "write_offset": None,
            }
        }

        self.lock_value = lock_value

        # Configure working states.
        self.init_methods_from_dict(
            {
                "import": {
                    "getter": None,
                    "processor": self.import_kmers,
                    "transition": self.finish_import,
                },
                "map": {
                    "getter": None,
                    "processor": self.map,
                    "transition": None,
                },
                "shutdown": self.save_shutdown
            }
        )

    @property
    def shm_keys(self):
        return self.shm_allocations["keys"]["array"]

    @shm_keys.setter
    def shm_keys(self, value):
        self.shm_allocations["keys"]["array"] = value

    @property
    def shm_values(self):
        return self.shm_allocations["values"]["array"]

    @shm_values.setter
    def shm_values(self, value):
        self.shm_allocations["values"]["array"] = value

    @property
    def shm_keys_allocation(self):
        return self.shm_allocations["keys"]["allocation"]

    @shm_keys_allocation.setter
    def shm_keys_allocation(self, value):
        self.shm_allocations["keys"]["allocation"] = value

    @property
    def shm_values_allocation(self):
        return self.shm_allocations["values"]["allocation"]

    @shm_values_allocation.setter
    def shm_values_allocation(self, value):
        self.shm_allocations["values"]["allocation"] = value

    @property
    def shm_keys_name(self):
        return self.shm_allocations["keys"]["shm_name"]

    @shm_keys_name.setter
    def shm_keys_name(self, value):
        self.shm_allocations["keys"]["shm_name"] = value

    @property
    def shm_values_name(self):
        return self.shm_allocations["values"]["shm_name"]

    @shm_values_name.setter
    def shm_values_name(self, value):
        self.shm_allocations["values"]["shm_name"] = value

    @property
    def shm_keys_size(self):
        return self.shm_allocations["keys"]["size"]

    @shm_keys_size.setter
    def shm_keys_size(self, value):
        self.shm_allocations["keys"]["size"] = value

    @property
    def shm_values_size(self):
        return self.shm_allocations["values"]["size"]

    @shm_values_size.setter
    def shm_values_size(self, value):
        self.shm_allocations["values"]["size"] = value

    @property
    def keys_data_type(self):
        return self.shm_allocations["keys"]["data_type"]

    @property
    def keys_data_type_size(self):
        return self.shm_allocations["keys"]["data_type_size"]

    @property
    def values_data_type(self):
        return self.shm_allocations["values"]["data_type"]

    @property
    def values_data_type_size(self):
        return self.shm_allocations["values"]["data_type_size"]

    @property
    def pile_length(self):
        return len(self.pile)

    @property
    def size(self):
        return sum([
            self.product(*kmers.shape, self.keys_data_type_size)
            for kmers in self.yield_kmers()
        ])

    @property
    def shape(self):
        return tuple(
            fn([kmers.shape[i] for kmers in self.yield_kmers()])
            for i, fn in enumerate((sum, max))
        )

    def yield_kmers(self):
        if self.shm_keys is not None:
            yield self.shm_keys

        for arr in self.pile:
            yield arr

    @staticmethod
    def is_int(string):
        try:
            value = float(string)
            int_value = int(value)

            if int_value != value:
                return False

            # Positive definite.
            if int_value < 1:
                return False

            return int_value

        except TypeError:
            return False

    @staticmethod
    def product(*args):

        total = 1
        for arg in args:
            total *= arg

        return total

    def import_kmers(self, data):
        window, shm_data = data
        kmers = self.view_kmers(window, shm_data)

        NO_KMERS_INFO = "No new kmers!"
        self.logger.info("Importing %d kmers at %d."
                         % (kmers.shape[0], self.pid))

        if kmers.shape[0] == 0:

            self.logger.info(NO_KMERS_INFO)

        else:

            # Make disjoint.
            kmers = self.disjoint(kmers)

            # If these kmers have no information to offer, continue.
            if kmers.shape[0] == 0:
                self.logger.info(NO_KMERS_INFO)

            else:
                self.pile.append(kmers)

                # Check if pile needs to be collapsed.
                if self.pile_length >= self.pile_threshold:
                    self.collapse()

                else:
                    self.logger.info("Process %d waiting this time around (%d/%s)."
                                     % (self.pid, self.pile_length, str(self.pile_threshold)))

    def view_kmers(self, window, data):
        start, end = window
        name, shape = data

        # Establish shared_memory connection.
        if name not in self.in_arrays:
            self.in_allocations[name] = multiprocessing.shared_memory.SharedMemory(name)
            self.in_arrays[name] = np.ndarray(
                shape=shape,
                dtype=self.dtypes.keys_dtype,
                buffer=self.in_allocations[name].buf
            )

        return (self.in_arrays[name])[start:end]

    @log_timer
    def disjoint(self, new_kmers):
        # Null case.
        n = new_kmers.shape[0]
        # Establish list of kmers against which to cull.
        comparisons = [self.shm_keys, *self.pile] if self.shm_keys is not None else self.pile

        # Remove any kmers that appear in the pile.
        n = make_disjoint([new_kmers, *comparisons])
        # Copy values into a new array.
        distinct_kmers = np.copy(new_kmers[:n])

        # Demonstrate completion.
        new_kmers[:] = 0

        return distinct_kmers

    def collapse(self):
        # Calculate new array parameters.
        size = self.size
        shape = self.shape

        # Remember old shape for filling.
        old_shape = 0 if self.shm_keys is None else self.shm_keys.shape[0]
        self.logger.info("Process %d sending extension request to shape (%d, %d)."
                         % (self.pid, shape[0], shape[1]))

        # Request extension to memory.
        create = self.shm_keys is None
        self.request_extension(size, create=create)

        # Remap memory extension.
        if self.shm_keys_allocation is None:
            self.shm_keys_allocation = ResizableSharedMemory(self.shm_keys_name)

        else:
            self.shm_keys_allocation.remap(size)

        self.shm_keys = np.ndarray(
            shape=shape,
            dtype=self.keys_data_type,
            buffer=self.shm_keys_allocation.buf
        )

        self.logger.info("Process %d extension successfully mapped." % self.pid)

        # Collapse the pile.
        self.collapse_pile(old_shape)

        self.logger.info("Process %d pile successfully collapsed." % self.pid)

    @log_timer
    def request_extension(self, size, create=False):
        self.con.send(
            (
                "shm_extend",
                (
                    self.shm_keys_name,
                    size,
                    create
                )
            )
        )
        result = self.con.recv()

        # Set shared_memory details.
        if create:
            self.shm_allocations["keys"]["shm_name"] = result
            self.shm_allocations["keys"]["active"] = True

        self.shm_keys_size = size

    @log_timer
    def collapse_pile(self, old_shape):
        size = self.pile_threshold + 1
        list_of_arrays = [arr for arr in self.yield_kmers()]

        # Check arrays are C contiguous.
        for arr in list_of_arrays:
            if not arr.flags["C_CONTIGUOUS"]:
                raise ValueError("Piled array is not C-contiguous!")

        filler(size, old_shape, list_of_arrays)

        # Reset the pile.
        self.pile = []

    def finish_import(self, *args):
        self.check_collapse()
        self.tell_size()

    def check_collapse(self):
        # Check if we need one final collapsing of the pile.
        if self.pile_length > 0:
            self.collapse()

        self.logger.info("Process %d finished." % self.pid)

    def tell_size(self):
        # Check for case where we have received no work.
        if self.shm_keys is not None and not self.is_supposed_to_die:
            self.con.send(
                (
                    "keys_size",
                    self.shape
                )
            )

            values_shm_name = self.con.recv()
            values_shape = self.shape[:1]

            # Initiate values allocation.
            self.shm_values_name = values_shm_name
            self.shm_values_size = self.product(values_shape, self.values_data_type_size)
            self.shm_values_allocation = ResizableSharedMemory(values_shm_name)
            self.shm_values = np.ndarray(
                shape=values_shape,
                dtype=self.values_data_type,
                buffer=self.shm_values_allocation.buf
            )

            self.shm_allocations["values"]["active"] = True

    def map(self, data):
        window, shm_data, node_id = data
        kmers = self.view_kmers(window, shm_data)

        self.logger.info("Process %d mapping %d." % (self.pid, node_id))
        binarysearch_put(
            kmers,
            self.shm_keys,
            self.shm_values,
            self.lca_matrix,
            node_id
        )

        # Demonstrate finished.
        kmers[:] = 0

    def save_shutdown(self, *args):
        self.save_data(*args)
        self.close_stores()
        self.close_shm()

    def save_data(self, *args):

        if self.shm_allocations["keys"]["active"] and self.shm_allocations["values"]["active"]:

            # Find our offset.
            self.con.send(
                (
                    "save",
                    self.shape
                )
            )

            arr_offset, db_path = self.con.recv()

            for arr in ("keys", "values"):
                self.shm_allocations[arr]["write_offset"] = arr_offset

            # Wait until it is my turn to open file.
            while self.loop_condition(self.lock_value.value != arr_offset):
                self.timeout()

            # Save data to disk.
            db = TablesDb(db_url=db_path, edit=True)

            for arr_name in ("keys", "values"):
                writer = db.write_keys if arr_name == "keys" else db.write_values
                arr = self.shm_allocations[arr_name]["array"]

                rows_written = writer(
                    arr=self.shm_allocations[arr_name]["array"],
                    offset=self.shm_allocations[arr_name]["write_offset"]
                )

                if rows_written != arr.shape[0]:
                    self.logger.info("Process %d writing of data failed!" % self.pid)
                else:
                    self.logger.info("Process %d successfully wrote data to %s!" % (self.pid, db_path))

            db.close()

            with self.lock_value.get_lock():
                self.lock_value.value += arr.shape[0]

        # Case where no data has been received.
        else:
            if not self.is_supposed_to_die:
                self.log("Warning! A child has gone this build without any work.")

                self.con.send(
                    (
                        "save",
                        None
                    )
                )

                receipt = self.con.recv()
                self.log(receipt)

    def close_stores(self):
        shm_checks = (
            (
                self.shm_allocations[col]["active"],
                self.shm_allocations[col]["allocation"]
            )
            for col in ("keys", "values")
        )

        for shm_check, shm_all in shm_checks:
            if shm_check:

                # Send an unlink and kill request to the main process.
                try:
                    shm_all.unlink()
                    shm_all.close()

                except FileNotFoundError:
                    pass

        # Log my death.
        self.logger.info("%d closed stores." % self.pid)

    def close_shm(self):
        for allocation in self.in_allocations.values():
            allocation.close()

