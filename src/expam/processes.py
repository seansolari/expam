import io
import math
import queue
import time
import multiprocessing.shared_memory
import traceback
from multiprocessing import Process, Queue, shared_memory
import os
import re

import numpy as np
import psutil
from expam.c.extract import get_raw_kmers, remove_duplicates, to_32bit, import_mask, map_kmers, \
    bytes_intersect
from expam.c.map import binarysearch, binarysearch_put, classify
from expam.c.sets import filler, disjoint as make_disjoint
from expam.run import log_timer, new_logger
from expam.run import CHUNK_SIZE, EXPAM_TEMP_EXT, KEYS_DATA_TYPE, KEYS_DATA_TYPE_SIZE, \
    MAP_DATA_TYPE, MAP_DATA_TYPE_SIZE, NULL_VALUE, TIMEOUT
from expam.sequences import format_name, get_opener
from expam.stores import MaskStore, ResizableSharedMemory, SequenceStore, TablesDb

COMMAND_CODES = {
    "TRANSITION": 2,
    "ACTIVE": 1,
    "SLEEP": 0,
    "DEAD": -1,
    "ERROR": -2
}


class JobWorker(Process):

    def __init__(self, logging_dir, status_index, statuses,
                 child_statuses, timeout=TIMEOUT, job_name="Worker",
                 phases=None):

        super(JobWorker, self).__init__()

        self._process = None

        self._status_id = status_index
        self._fellow_statuses = statuses
        self._child_statuses = child_statuses

        self.phases = [] if phases is None else list(phases)
        self.work_queue = None  # Where work comes from.
        self.child_queues = None  # Where to distribute work.

        #
        # Defined at ControlCenter.
        self._parents = []
        self._children = []
        self._command_permissions = {
            "command": "None",
            "counts": 0
        }

        #
        # User defined methods.
        self._getters = {}
        self._processors = {}
        self._transitions = {}
        self._shutdown = self.default_shutdown

        self._timeout = timeout
        self.job_name = job_name

        # Logging to disk.
        self.logging_dir = logging_dir
        self.logger = None

    def set_method(self, method_name, phase, method):
        methods = {
            "getter": self._getters,
            "processor": self._processors,
            "transition": self._transitions,
        }

        if method_name not in methods:
            raise ValueError("Invalid method name %s!" % method_name)

        if method is not None:
            methods[method_name][phase] = method

    def init_methods_from_dict(self, config):

        if not isinstance(config, dict):
            raise TypeError("Config must be dictionary not %s!" % str(type(config)))

        for key, value in config.items():
            if value is None:
                continue

            if key == "shutdown":

                if not callable(value):
                    raise TypeError("Shutdown object must be callable!")

                self.set_shutdown(value)

            elif key in self.phases:

                if not isinstance(value, dict):
                    raise TypeError("Method declaration must be dictionary not %s!"
                                    % str(type(value)))

                for method_name, method in value.items():
                    self.set_method(
                        method_name=method_name,
                        phase=key,
                        method=method
                    )

            else:
                raise ValueError("Unknown key in class config: %s!" % key)

    def set_shutdown(self, method):
        self._shutdown = method

    def children_have_obeyed(self, *codes):
        if self._child_statuses is None:
            return True

        return all([
            child_status in codes
            for child_status in self._child_statuses
        ])

    @property
    def is_primary_child(self):
        return True if self._status_id == 0 else False

    @property
    def children_have_errors(self):
        if self._child_statuses is None:
            return False

        return any([
            child_status == COMMAND_CODES["ERROR"]
            for child_status in self._child_statuses
        ])

    @property
    def error_propagated(self):
        if self._child_statuses is None:
            return True

        return all([
            child_status == COMMAND_CODES["ERROR"]
            for child_status in self._child_statuses
        ])

    @property
    def children_dead(self):
        return self.children_have_obeyed(COMMAND_CODES["DEAD"])

    def peers_are_okay(self, cid):
        return all([
            peer_status != COMMAND_CODES["ERROR"]
            for i, peer_status in enumerate(self._fellow_statuses)
            if i != cid
        ])

    def wake_children(self):
        for cid in self.yield_children_id():
            self._child_statuses[cid] = COMMAND_CODES["ACTIVE"]

    @property
    def child_id(self):
        return self._status_id

    @property
    def n_parents(self):
        return len(self._parents)

    @property
    def n_children(self):
        return len(self._children)

    @property
    def n_peers(self):
        return len(self._fellow_statuses)

    @property
    def phase(self):
        return self.phases[0]

    @property
    def status(self):
        return self._fellow_statuses[self._status_id]

    @property
    def name(self):
        return self.job_name

    @status.setter
    def status(self, value):
        assert value in COMMAND_CODES.values()

        self._fellow_statuses[self._status_id] = value

    @property
    def current_command(self):
        return self._command_permissions["command"]

    @current_command.setter
    def current_command(self, value):
        self._command_permissions["command"] = value

    @property
    def command_permissions(self):
        return self._command_permissions["counts"]

    @command_permissions.setter
    def command_permissions(self, value):
        self._command_permissions["counts"] = value

    @property
    def is_supposed_to_die(self):
        return self.status == COMMAND_CODES["ERROR"]

    @property
    def is_supposed_to_sleep(self):
        return self.status == COMMAND_CODES["SLEEP"]

    def reset_command_permissions(self):
        self.current_command = "None"
        self.command_permissions = 0

    def tell_children(self, *args):
        if len(args) == 1:
            args = args * self.n_children

        assert len(args) == self.n_children

        for i in range(self.n_children):
            self.tell_child(i, args[i])

    def tell_child(self, cid, data):

        if not self.child_is_dead(cid):
            self.child_queues[cid].put(data)

    def message_child(self, cid, err_data):
        self._children[cid].send(err_data)

    def yield_children_id(self):
        for i in range(len(self._children)):
            yield i

    def timeout(self):
        time.sleep(self._timeout)

    def is_active(self):
        if self.status == COMMAND_CODES["ACTIVE"]:
            return True
        return False

    def _increment_command(self, command_string):

        if self.current_command == "None":
            self.current_command = command_string

        elif self.current_command != command_string:
            raise ValueError("New command %s interrupted old command flow (%s)!"
                             % (command_string, self.current_command))

        # Check completion condition.
        self.command_permissions += 1

        if self.command_permissions == self.n_parents:

            command_code = COMMAND_CODES[command_string]
            death_code = COMMAND_CODES["DEAD"]

            # Broadcast sleep command to all children.
            self.log("%d moving command to children..." % self.pid)
            self.tell_children((command_code, command_string))
            self.log("%d children told!" % self.pid)

            # Children must have either slept or have died due to no more phases
            while not self.children_have_obeyed(command_code, death_code):
                if self.is_primary_child:
                    has_errors = self.check_children()  # Monitor for errors.

                    if has_errors:
                        return

                # Listen for messages.
                messages = self.check_for_messages()
                if messages:
                    self.begin_error_sequence(messages)

                    return

                self.timeout()

            self.log("%d my children have obeyed!" % self.pid)

            # Start transition to sleep state.
            self.reset_command_permissions()
            self.status = COMMAND_CODES["TRANSITION"]

    def activate(self):
        if not self.is_active():
            self.status = COMMAND_CODES["ACTIVE"]

    def process(self, *args):
        memory_info = self._process.memory_info()
        self.log("%d currently using (RSS: %fGb, VSZ: %fGb) of RAM. Processing data..."
                 % (self.pid,
                    memory_info.rss / 1024 ** 3,
                    memory_info.vms / 1024 ** 3))

        return self._processors[self.phase](*args)

    def sleep(self, *args):
        self.log("%d received sleep command." % self.pid)

        self._increment_command("SLEEP")

        if self.status == COMMAND_CODES["TRANSITION"]:
            self.transition(*args)

            # Sleep when phase transition is complete.
            self.status = COMMAND_CODES["SLEEP"]

            # Wait to be awoken by parent.
            while not self.is_active():
                self.timeout()

            # Wake up children.
            self.wake_children()

    def die(self, *args):
        self.log("%d received die command." % self.pid)

        self._increment_command("DEAD")

        if self.status == COMMAND_CODES["TRANSITION"]:
            self.shutdown(*args)

            # Die when phase transition is complete.
            self.status = COMMAND_CODES["DEAD"]

    def set_error_status(self):
        self.status = COMMAND_CODES["ERROR"]

    def reset_status(self):
        self.status = COMMAND_CODES["ACTIVE"]

    def set_dead_status(self):
        self.status = COMMAND_CODES["DEAD"]

    def die_error(self, *args):
        self.set_error_status()

        # Clear my work queue to shut down safely.
        self.clear_work_queue()

        self.shutdown()

    def clear_work_queue(self):
        while True:
            try:
                _ = self.work_queue.get(block=True, timeout=self._timeout)

            except queue.Empty:
                return

    def get_work(self, *args):

        # First check for any errors raised.
        msg = self.check_for_messages()

        if msg is not False:  # This is an error to be propagated.
            return msg

        else:
            if self.phases and self.phase in self._getters:
                return self._getters[self.phase](*args)

            assert len(args) == 2
            wait, pid = args[0], args[1]

            # Check for more jobs.
            try:
                work = self.work_queue.get(block=True, timeout=self._timeout)

                return work

            except queue.Empty:
                if not wait:
                    return None

    def loop_condition(self, cond):
        if self.received_message():
            return False

        return cond

    def received_message(self):
        if self._parents[0].poll(self._timeout):
            return True

        return False

    def check_for_messages(self):
        if self._parents[0].poll(self._timeout):
            return self._parents[0].recv()

        else:
            return False

    def transition(self, *args):
        if self.phases and self.phase in self._transitions:
            self._transitions[self.phase](*args)

        self.phases.pop(0)
        self.log("%d now at phase %s!" % (self.pid, ("end" if not self.phases else self.phase)))

    def shutdown(self, *args):
        return self._shutdown(*args)

    def default_shutdown(self, *args):
        self.log("%d %s closing..." % (self.pid, self.job_name))

    def raise_error(self, message):
        # If we already have an error message, that takes precedence.
        # Only send if there are no current messages (ie. peers haven't raised).

        self.set_error_status()

        if self._parents[0].poll(self._timeout):
            self.log("%d error was overridden!" % self.pid)

            # Another error was raised.
            self.begin_error_sequence()

        else:
            self._parents[0].send(message)
            self.log("%d error raised..." % self.pid)

            # Await confirmation that my error was received.
            self.log("%d awaiting response..." % self.pid)
            receipt = self._parents[0].recv()
            self.log("%d response received!" % self.pid)

            if receipt is True:
                self.log("%d sent error message!" % self.pid)

                if self.is_primary_child:
                    self.propagate_error(-1, message)

                self.die_error(message[1])

            else:
                self.log("%d error was overridden!" % self.pid)

                # Another error was raised first.
                self.begin_error_sequence(receipt)

    def propagate_error(self, eid, err_data):
        self.log("%d propagating error..." % self.pid)

        # Tell children about error and wait for them to shut down.
        for cid in self.yield_children_id():
            if cid != eid:
                self.message_child(cid, err_data)

        while not self.children_dead:
            self.timeout()

        self.log("%d error propagated!" % self.pid)

    def check_children(self):
        assert self.is_primary_child

        if self.children_have_errors:
            self.log("%d found errors!" % self.pid)
            self.set_error_status()

            # Get error and data.
            for cid in self.yield_children_id():

                if self._children[cid].poll(self._timeout):
                    err_data = self._children[cid].recv()
                    self._children[cid].send(True)  # Send confirmation.

                    break

            else:
                raise Exception("No errors found even though children are dead!")

            self.propagate_error(cid, err_data)

            # Tell parent.
            self.raise_error(err_data)

            return True

        return False

    def child_is_dead(self, cid):
        if self._child_statuses[cid] == COMMAND_CODES["DEAD"]:
            return True

        return False

    def begin_error_sequence(self, err_data=None):
        if err_data is None:
            err_data = self._parents[0].recv()

        if self.is_primary_child:
            self.propagate_error(-1, err_data)

        self.die_error(err_data)

    def log(self, msg):
        self.logger.info(msg)

    def run(self):
        # Activate logging.
        self.logger = new_logger(self.logging_dir, str(self.pid) + "_" + self.job_name)
        self.log("%d %s is now awake!" % (self.pid, self.job_name))

        # Activate memory tracking.
        self._process = psutil.Process(self.pid)

        # Phases loop.
        while not self.is_supposed_to_die and self.phases:

            # Check for errors.
            if self.is_primary_child:
                has_errors = self.check_children()

                if has_errors:
                    break

            # Processing loop.
            for pid in range(self.n_parents):

                try:
                    next_job = self.get_work(False, pid)  # True ==> wait for work.

                    if next_job is None:
                        break

                    command_code, data = next_job

                    if command_code == COMMAND_CODES["ACTIVE"]:
                        self.process(data)

                    elif command_code == COMMAND_CODES["SLEEP"]:
                        self.sleep(data)

                        if self.is_supposed_to_sleep:
                            break  # Check next phase exists.

                    elif command_code == COMMAND_CODES["DEAD"]:
                        self.die(data)

                        if self.is_supposed_to_die:
                            break  # Check next phase exists.

                    elif command_code == COMMAND_CODES["ERROR"]:
                        # Error message is currently propagating.
                        self.begin_error_sequence((command_code, data))

                        break

                    else:
                        raise ValueError("Invalid command code %d!" % command_code)

                except Exception as e:
                    # Log error with context.
                    self.log("*** Exception raised! ***\n" + traceback.format_exc())

                    self.log("Exception encountered in %d! Attempting raise..." % self.pid)

                    err_data = (COMMAND_CODES["ERROR"], e)
                    self.raise_error(err_data)

                    break

        if not self.is_supposed_to_die:
            self.shutdown()

        self.set_dead_status()


class UnionWorker(JobWorker):
    def __init__(self, main_con, lca_matrix, pile_size, logging_dir,
                 lock_value, status_index, statuses, child_statuses,
                 timeout=TIMEOUT, phases=None):

        super(UnionWorker, self).__init__(logging_dir, status_index, statuses, child_statuses,
                                          timeout=timeout, phases=phases)

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
        self.shm_allocations = {
            "keys": {
                "array": None,
                "allocation": None,
                "shm_name": None,
                "size": None,
                "data_type": KEYS_DATA_TYPE,
                "data_type_size": KEYS_DATA_TYPE_SIZE,
                "active": False,
                "write_offset": None,
            },
            "values": {
                "array": None,
                "allocation": None,
                "shm_name": None,
                "size": None,
                "data_type": MAP_DATA_TYPE,
                "data_type_size": MAP_DATA_TYPE_SIZE,
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
                dtype=KEYS_DATA_TYPE,
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


class ExtractWorker(JobWorker):
    def __init__(self, k, shm_kmer_params, kmer_ranges, node_map, out_dir, logging_dir,
                 status_index, statuses, child_statuses, timeout=TIMEOUT,
                 phases=None):

        super(ExtractWorker, self).__init__(logging_dir, status_index, statuses, child_statuses,
                                            timeout=timeout, phases=phases)

        self.job_name = "extract"
        self.k = k
        self.kmer_ranges = kmer_ranges

        self.chunks = (k // CHUNK_SIZE) + (k % CHUNK_SIZE != 0)
        self.d_type = np.dtype(
            [
                ('index', KEYS_DATA_TYPE),
                *(
                    ('k%d' % chunk, KEYS_DATA_TYPE)
                    for chunk in range(self.chunks)
                )
            ]
        )

        # For mapping LCAs.
        self.node_map = node_map

        # Hold the sequences and masks that I have extracted.
        self.sequence_store = SequenceStore(k, out_dir)
        self.mask_store = MaskStore()

        shm_kmer_name, shm_kmer_shape = shm_kmer_params
        allocation_shape = (shm_kmer_shape[0], shm_kmer_shape[1] + 1)

        # Long-term arrays used for extraction.
        self.allocation = np.zeros(
            shape=allocation_shape,  # Account for index col.
            dtype=KEYS_DATA_TYPE,
            order="C",  # Row-major storage (consecutive k-mer access).
        )

        # Long-term array for sending kmers to children.
        self.out_shape = shm_kmer_shape
        self.out_name = shm_kmer_name
        self.out_allocation = multiprocessing.shared_memory.SharedMemory(shm_kmer_name)
        self.out_array = np.ndarray(
            shape=shm_kmer_shape,
            dtype=KEYS_DATA_TYPE,
            buffer=self.out_allocation.buf
        )

        # Primary queue for getting work during mapping phase.
        self.map_queue = Queue()

        # Initialise methods for execution.
        self.init_methods_from_dict(
            {
                "import": {
                    "getter": None,
                    "processor": self.import_sequence,
                    "transition": self.prepare_map,
                },
                "map": {
                    "getter": self.get_my_sequence,
                    "processor": self.send_for_mapping,
                    "transition": self.empty_queue,
                },
                "shutdown": self.close_shm
            }
        )

    @log_timer
    def import_sequence(self, sequence_url):
        def resize_base(size):
            nonlocal base

            base = base[:size]
            update_view()

        def new_view():
            nonlocal base
            return base.ravel().view(dtype=self.d_type)

        def update_view():
            nonlocal base, view
            view = base.ravel().view(dtype=self.d_type)

        # Collect raw sequence.
        self.log("Importing sequence at %s." % sequence_url)
        seq_id, seq_length = self.sequence_store.add_sequence(sequence_url)

        # Get view on required segment of allocation.
        base = self.allocation[:seq_length]
        view = new_view()

        # Create draft mask and get corresponding kmers.
        new_size = get_raw_kmers(
            sequence=self.sequence_store[seq_id],
            k=self.k,
            arr=base
        )
        resize_base(new_size)

        # Sort by kmers, which also sorts the indices
        view.sort(
            order=['k%d' % chunk for chunk in range(self.chunks)],
            kind="stable"
        )
        self.log("Removing duplicates.")
        new_size = remove_duplicates(base)
        resize_base(new_size)

        # Save mask as 32-bit array.
        self.mask_store[seq_id] = to_32bit(view['index'])
        self.send_kmers(
            ranges=self.kmer_ranges,
            kmers=base[:, 1:]
        )

    @staticmethod
    def clear_temp_arrays(arr):
        arr[:] = 0

    @log_timer
    def send_kmers(self, ranges, kmers, node_id=None):
        target = np.ndarray(shape=(kmers.shape[1]), dtype=kmers.dtype)
        target[:] = ranges[0][0]

        total = 0
        n = len(ranges) - 1
        lower_ind = upper_ind = self.get_index(target, kmers, greater=False)

        # Wait for previous data to be sent.
        while self.loop_condition(not np.all(self.out_array == 0)):
            self.timeout()

        # Place kmers in out_array.
        self.out_array[:kmers.shape[0], :] = kmers[:]

        # Separate data between children.
        for i, (lower_bound, upper_bound) in ranges.items():
            target[:] = upper_bound
            upper_ind = self.get_index(target, kmers)

            # Upper bound can only attain equality for maximum hash value.
            if i == n:
                while (upper_ind < kmers.shape[0]) and self.kmers_equal(target, kmers[upper_ind]):
                    upper_ind += 1

            # Don't send empty work.
            if lower_ind != upper_ind:

                message_data = (
                    (lower_ind, upper_ind),
                    (self.out_name, self.out_shape),
                )
                if node_id is not None:
                    message_data += (node_id,)

                self.tell_child(
                    cid=i,
                    data=(
                        COMMAND_CODES["ACTIVE"],
                        message_data
                    )
                )

                total += upper_ind - lower_ind

            # Beginning next search.
            lower_ind = upper_ind

        self.log("%d / %d kmers sent!" % (total, kmers.shape[0]))

    def get_index(self, target, kmers, greater=True):
        def _check_ind():
            nonlocal ind, kmers

            if (greater and ind == kmers.shape[0]) or ((not greater) and ind == 0):
                return True

            return False

        # Find closest element in list. Note that this element
        # may be less than the bound.
        _, ind = binarysearch(target, kmers)

        if _check_ind():
            return int(ind)

        while (not _check_ind()) and self.kmers_compare(kmers[ind], target, greater):
            if greater:
                ind += 1

            else:
                ind -= 1

        return int(ind)

    @staticmethod
    def kmers_equal(one, two):
        return np.all(one == two)

    @staticmethod
    def kmers_compare(one, two, greater):
        for i in range(one.shape[0]):
            if one[i] < two[i]:
                return greater

            elif one[i] > two[i]:
                return not greater

            elif i + 1 == one.shape[0]:
                return not greater

    def prepare_map(self, *args):
        # Save accession ids from my sequences.
        self.sequence_store.save_accession_ids(suffix=str(self._status_id))

        # Prepare mapping jobs.
        self.load_queue()

    def load_queue(self):
        for seq_id in self.sequence_store.keys():
            self.map_queue.put(
                (
                    COMMAND_CODES["ACTIVE"],
                    seq_id
                )
            )

    def get_my_sequence(self, *args):

        # Try to get local work first.
        try:
            work = self.map_queue.get(block=True, timeout=self._timeout)

            return work

        except queue.Empty:
            pass

        # Now look for commands from parent.
        try:
            work = self.work_queue.get(block=True, timeout=self._timeout)

            return work

        except queue.Empty:
            return None

    def send_for_mapping(self, seq_id):
        node_id = self.node_map[seq_id]

        sequence = self.sequence_store[seq_id]
        mask = self.mask_store[seq_id]
        n_kmers = mask.shape[0]

        # Get view on required segment of allocation.
        base = self.allocation[:n_kmers]
        view = base.ravel().view(dtype=np.dtype(
            [
                ('index', KEYS_DATA_TYPE),
                *(
                    ('k%d' % chunk, KEYS_DATA_TYPE)
                    for chunk in range(self.chunks)
                )
            ]
        ))

        # Place 32-bit mask in 64-bit array.
        import_mask(
            data=mask,
            dest=view['index']
        )
        map_kmers(
            sequence=sequence,
            k=self.k,
            arr=base
        )

        # Send these kmers to be mapped.
        self.send_kmers(
            ranges=self.kmer_ranges,
            kmers=base[:, 1:],
            node_id=node_id
        )

    def empty_queue(self, *args):
        try:
            _ = self.map_queue.get(block=True, timeout=self._timeout)

            raise ValueError("Mapping queue was not empty!")

        except queue.Empty:

            return

    def close_shm(self):
        try:
            self.out_allocation.unlink()
            self.out_allocation.close()

        except FileNotFoundError:
            pass


class ReadWorker(JobWorker):
    def __init__(self, paired_mode, n_classifiers, logging_dir, status_index, statuses, child_statuses,
                 timeout=TIMEOUT, phases=None):
        super(ReadWorker, self).__init__(logging_dir, status_index, statuses, child_statuses, timeout=timeout,
                                         phases=phases)
        self.job_name = "reader"

        self.paired_mode = paired_mode

        self.n_classifiers = n_classifiers
        self._c = 0  # Which classifier do we next send data to?

        self.n_streams = 1 if not paired_mode else 2

        self.init_methods_from_dict(
            {
                "classify": {
                    "getter": None,
                    "processor": self.load_reads,
                    "transition": None,
                }
            }
        )

    def load_reads(self, file_urls):
        print("Loading reads from %s..." % ", ".join(file_urls))

        file_names, file_modes, file_openers, file_types = self.init_openers(*file_urls)

        # Open each file.
        _file_handles = [
            opener(file_url, mode)
            for opener, file_url, mode in zip(file_openers, file_urls, file_modes)
        ]

        read_stream = io.BytesIO()

        read_ids = []
        read_seq = None

        while read_ids or read_seq is None:  # read_seq is None is initial case.
            for i in range(self.n_streams):
                file_type = file_types[i]

                if file_type == "fasta":
                    next_read_id = self.scan_fa_to_next_read(_file_handles[i], read_stream)

                else:  # File type is "fastq".
                    next_read_id = self.scan_fq_to_next_read(_file_handles[i], read_stream)

                if next_read_id is not None:
                    if i == 0:
                        read_ids.append(next_read_id)
                    else:
                        read_ids[-1] = bytes_intersect(read_ids[-1], next_read_id)

            read_seq = read_stream.getvalue()

            if len(read_seq) > self.n_streams:
                read_id = read_ids.pop(0)

                self.send_read(min(file_names), read_id, read_seq)

            read_stream.seek(0)
            read_stream.truncate()

        for i in range(self.n_streams):
            # Close file handles.
            _file_handles[i].close()

        # Close data stream.
        del read_stream

    def send_read(self, file_name, read_id, read):
        self.tell_child(
            cid=self._c,
            data=(
                COMMAND_CODES["ACTIVE"],
                (file_name, read_id, read, self.child_id)
            )
        )

        self._c = (self._c + 1) % self.n_classifiers

    def init_openers(self, *file_urls):
        file_names, file_modes, file_openers, file_types = [], [], [], []

        for file_url in file_urls:
            file_name = format_name(file_url)
            file_name, (mode, opener) = get_opener(file_name)

            if self.check_suffix(file_url, ".fasta", ".fa"):
                file_type = "fasta"

            elif self.check_suffix(file_url, ".fastq", ".fq"):
                file_type = "fastq"

            else:
                raise ValueError("Unknown sequence format for %s!" % file_url)

            file_names.append(file_name)
            file_modes.append(mode)
            file_openers.append(opener)
            file_types.append(file_type)

        return file_names, file_modes, file_openers, file_types

    @staticmethod
    def check_suffix(string, *suffixes):
        for suffix in suffixes:
            if re.findall(r"(%s)(?:.tar)*(?:.gz)*$" % suffix, string):
                return True

        return False

    @staticmethod
    def scan_fa_to_next_read(file, read_stream):
        next_read_id = None

        while True:
            line = file.readline()

            if line == b'':
                read_stream.write(b'N')
                break  # EOF.

            if b'>' in line:
                next_read_id = re.match(b'^>(\S+)\s', line).group(1)

                # End my contribution to this read.
                read_stream.write(b'N')
                break

            else:
                read_stream.write(line.strip())

        return next_read_id

    @staticmethod
    def scan_fq_to_next_read(file, read_stream):
        next_read_id = None

        while True:
            line = file.readline()

            if line == b'':
                read_stream.write(b'N')
                break

            if b'@' in line:
                next_read_id = re.match(b'^@(.*)\s', line).group(1)

                read_stream.write(b'N')
                break

            else:
                read_stream.write(line.strip())

                # Skip next two lines.
                for _ in range(2):
                    file.readline()

        return next_read_id


class ClassifyWorker(JobWorker):
    def __init__(self, k, out_dir, temp_dir, kmer_db, lca_matrix, n_readers, logging_dir, alpha, status_index,
                 statuses, child_statuses, timeout=TIMEOUT, phases=None):
        super(ClassifyWorker, self).__init__(logging_dir, status_index, statuses, child_statuses, timeout=timeout,
                                             phases=phases)
        self.job_name = "classify"
        self.n_readers = n_readers

        # Configuration parameters.
        self.k = k
        self.alpha = alpha
        self.chunks = (k // CHUNK_SIZE) + (k % CHUNK_SIZE != 0)

        self.out_dir = out_dir
        self.temp_dir = temp_dir

        self._read_array_length = 150
        self.read_array = np.ndarray((self._read_array_length, self.chunks + 1), dtype=KEYS_DATA_TYPE)

        self.RAW = "_%s.reads_%s" % ("%d", EXPAM_TEMP_EXT)

        # Load the database.
        self.keys_shm = shared_memory.SharedMemory(
            name=kmer_db.keys_shm_name,
            create=False
        )
        self.keys = np.ndarray(
            shape=kmer_db.keys_meta["shape"],
            dtype=kmer_db.keys_meta["dtype"],
            buffer=self.keys_shm.buf
        )
        self.values_shm = shared_memory.SharedMemory(
            name=kmer_db.values_shm_name,
            create=False
        )
        self.values = np.ndarray(
            shape=kmer_db.values_meta["shape"],
            dtype=kmer_db.values_meta["dtype"],
            buffer=self.values_shm.buf
        )

        self.current_file_data = [("", []) for _ in range(n_readers)]

        # LCA Processing.
        self.lca_matrix = lca_matrix

        self.init_methods_from_dict(
            {
                "classify": {
                    "getter": None,
                    "processor": self.classify_read,
                    "transition": self.save_final_data,
                }
            }
        )

    def classify_read(self, read_data):
        new_file_name, read_id, read, reader_id = read_data

        # Continue to get reads without checking for errors; only check for errors when there is a break in workflow.
        while True:
            read_id = read_id.decode("utf8")

            max_n_kmers = len(read) - self.k + 1
            if max_n_kmers > self.read_array.shape[0]:
                self.read_array = np.ndarray((max_n_kmers, self.chunks + 1), dtype=KEYS_DATA_TYPE)

            if new_file_name != self.current_file_data[reader_id][0]:
                if len(self.current_file_data[reader_id][1]) > 0:
                    self.save_temp_data(reader_id)

                self.current_file_data[reader_id] = (new_file_name, [])

            result, best_clade, read_length, result_string = self.process_read(read)
            self.current_file_data[reader_id][1].append(
                [result, read_id, "p" + str(best_clade), str(read_length), result_string]
            )

            # Get a head start on any available work.
            try:
                command_code, data = self.work_queue.get(block=True, timeout=self._timeout)

                if command_code != COMMAND_CODES["ACTIVE"]:
                    self.work_queue.put((command_code, data))
                    break

                new_file_name, read_id, read, reader_id = data
            except queue.Empty:
                break

    def save_temp_data(self, reader_id):
        def format_result(result_list):
            return "\n".join(["\t".join(item) for item in result_list])

        def base_name(file_url):
            return re.match(r"(\S+)\.(?:fa|fq|fasta|fastq)(?:\.tar)*(?:\.gz)*$", file_url).group(1)

        file_name, file_data = self.current_file_data[reader_id]

        my_file_name = os.path.join(self.temp_dir, base_name(file_name) + (self.RAW % self.pid))
        with open(my_file_name, "a+") as f:
            f.write(format_result(file_data))

    def process_read(self, read):
        read_length = len(read)
        n_kmers = get_raw_kmers(read, self.k, self.read_array)

        # Control for reads that don't have long enough good regions for kmers to be extracted.
        if n_kmers == 0:
            return "U", 0, read_length, ""

        # Map the kmers.
        cls, common_clade, kmer_string = classify(self.read_array[:n_kmers], self.keys, self.values,
                                                  self.lca_matrix, NULL_VALUE, self.alpha)

        if cls == -1:
            result = "U"

        elif cls == 1:
            result = "S"

        else:
            result = "C"

        return result, common_clade, read_length, kmer_string

    def save_final_data(self, *args):
        for reader_id in range(self.n_readers):
            if len(self.current_file_data[reader_id][1]) > 0:
                self.save_temp_data(reader_id)

        self.keys_shm.close()
        self.values_shm.close()
