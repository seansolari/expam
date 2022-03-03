import datetime
import json
import logging
import subprocess

import math
import os
import queue
import time
import sys
from ctypes import *
from multiprocessing import Array, Pipe, Queue, shared_memory, Value
from timeit import default_timer

import matplotlib.pyplot as plt
import expam
import numpy as np

KEYS_DATA_TYPE = np.uint64
KEYS_DATA_TYPE_SIZE = np.dtype(KEYS_DATA_TYPE).itemsize
KEYS_DTYPE_STR = "uint64"

MAP_DATA_TYPE = np.uint16
MAP_DATA_TYPE_SIZE = np.dtype(MAP_DATA_TYPE).itemsize
MAP_DTYPE_STR = "uint16"

CHUNK_SIZE = 32
UNION_RATIO = 4 / 5

TIMEOUT = 1e-4

ACCESSION_ID_PATH = os.path.join("phylogeny", "accession_ids.csv")
COORDS_PATH = os.path.join("phylogeny", "coords.json")
LCA_MATRIX_PATH = os.path.join("phylogeny", "lca_matrix.npy")
DB_PATH = os.path.join("database", "expam_db.h5")

NULL_VALUE = 0  # Return value for foreign kmer.
EXPAM_TEMP_EXT = "expam_temp"

LOG_FORMAT = logging.Formatter(
    fmt='%(asctime)s... %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p'
)


def die(msg):
    print("ERROR\n", msg, "\n")
    sys.exit(1)


class Timer:
    def __init__(self):
        self.start = None
        self.elapsed = None

    def __enter__(self):
        self.start = default_timer()
        return self

    def __exit__(self, exit_type, value, traceback):
        self.elapsed = default_timer() - self.start

    def __str__(self):
        return self.verbose()

    def __float__(self):
        return self.elapsed

    def verbose(self):
        if self.elapsed is None:
            return "<not_measured>"
        return str(self.elapsed) + "s"


class JSONConfig:
    _parse_data_type = {
        "uint8": np.uint8,
        "uint16": np.uint16,
        "uint32": np.uint32,
        "uint64": np.uint64,
        "int8": np.int8,
        "int16": np.int16,
        "int32": np.int32,
        "int64": np.int64,
    }

    def __init__(self, url=None, name=""):
        self.url = url

        # Default values.
        self.params = {
            "name": name,
            "groups": {
                'default': {
                    'sequence_files': [],
                    'phylogeny_path': [],
                    'k': None,
                    's': None,
                },
            },
            "phylogeny_path": None,
            "k": None,
            "n": None,
            "s": None,
            "pile": None,
            "keys_shape": None,
            "keys_data_type": None,
            "values_shape": None,
            "values_data_type": None,
        }

        if url is not None:
            with open(url, 'r') as f:
                _conf = json.load(f)

                if not isinstance(_conf, dict):
                    die('Invalid configuration file!')

                for key in _conf:
                    if key == 'groups':
                        self.params['groups'].update(_conf['groups'])

                    else:
                        self.params[key] = _conf[key]

    def __getitem__(self, key):
        return self.params[key]

    def __str__(self):
        name_line = "<<< expam configuration file: %s >>>\n" % self.params['name']
        data_line = "phylogeny\t-->\t%s\nk\t\t-->\t%s\nn\t\t-->\t%s" % (self.params['phylogeny_path'],
                                                                    self.params['k'], self.params['n'])
        par_line = "sketch\t\t-->\t%s\npile\t\t-->\t%s" % (self.params['s'], self.params['pile'])

        group_fmt = "\n----------------\ngroup name: %s\n\tk\t\t-->\t%s\n\tsketch\t\t-->\t%s\n\tsequences\t-->\t%d\n"
        group_line = "\n".join([
            group_fmt % (name, group['k'], group['s'], len(group['sequence_files']))
            for name, group in self.params['groups'].items()
        ])

        return "\n".join([name_line, data_line, par_line, group_line])

    def get_paths(self):
        return [
            genome_path
            for group in self.params["groups"].values()
            for genome_path in group['sequence_files']
        ]

    def new_group(self, name, files=None, phylogeny_path=None, k=None, s=None):
        if name in self.params['groups']:
            raise ValueError("Group %s already exists!" % name)

        self.params['groups'][name] = {
            'sequence_files': [] if files is None else files,
            'phylogeny_path': phylogeny_path,
            'k': k,
            's': s,
        }

    def set(self, **kwargs):
        for arg in kwargs:
            if arg in self.params and arg != 'groups':
                if kwargs[arg] is not None:
                    self.params[arg] = kwargs[arg]

            else:
                raise ValueError("Unknown configuration category %s!" % arg)

    def group_set(self, group_name, **kwargs):
        for arg in kwargs:
            if arg in self.params['groups'][group_name] and arg != 'sequence_files' and kwargs[arg] is not None:
                self.params['groups'][group_name][arg] = kwargs[arg]

    def group_get(self, group_name):
        if group_name not in self.params['groups']:
            raise IndexError("Invalid group %s!" % group_name)

        group = self.params['groups'][group_name]
        attrs = {attr: group[attr] for attr in ('k', 's', 'sequence_files')}

        for req_param in ('k', 's'):
            if group[req_param] is None:
                if self.params[req_param]:
                    attrs[req_param] = int(self.params[req_param])

            else:
                attrs[req_param] = int(attrs[req_param])

        return attrs['k'], attrs['s'], attrs['sequence_files']

    def add_sequence(self, file_name, group=None):
        if group is None:
            group = 'default'

        if group not in self.params['groups']:
            self.new_group(group)

        if not isinstance(file_name, list):
            file_name = {file_name}
        else:
            file_name = set(file_name)

        # Check these sequences aren't already added.
        for current_group in self.params['groups']:
            if current_group != group:
                current_seqs = set(self.params['groups'][current_group]['sequence_files'])

                if file_name & current_seqs:
                    die("Some of these sequences are already in the database!")

        self.params['groups'][group]['sequence_files'].extend(file_name)

    def remove_sequence(self, file_name, group=None):
        if group is None:
            group = 'default'

        if isinstance(file_name, list):
            file_name = set(file_name)
        else:
            file_name = {file_name}

        sequences = set(self.params['groups'][group]['sequence_files'])
        sequences -= file_name

        self.params['groups'][group]['sequence_files'] = list(sequences)

    def remove_group(self, group=None):
        if group is None:
            self.params['groups']['default'] = {
                    'sequence_files': [],
                    'phylogeny_path': [],
                    'k': None,
                    's': None,
                }
        elif group in self.params['groups']:
            self.params['groups'].pop(group)
        else:
            raise Exception("No such group %s!" % group)

    def _dump(self):
        return json.dumps(self.params)

    def groups(self):
        return list(self.params['groups'].keys())

    def save(self, url=None):
        if url is None:
            if self.url is None:
                raise Exception("Save url required!")

            url = self.url

        with open(url, 'w') as f:
            json.dump(self.params, f)


class ControlCenter:
    def __init__(self, logging_dir, group_name=None, workers=None, child_statuses=None, phases=None,
                 phase_queues=None, child_queues=None, children=None, processors=None, transitions=None,
                 timeout=TIMEOUT):

        self.group_name = group_name
        self.workers = [] if workers is None else list(workers)

        self.error = False  # Error state.
        self._child_statuses = child_statuses  # Shared array.

        self.phases = [] if phases is None else list(phases)
        self.phase_queues = {} if phase_queues is None else dict(phase_queues)
        self.child_queues = child_queues

        # Connection objects.
        self._children = [] if children is None else list(children)

        # User defined methods.
        self._processors = {} if processors is None else dict(processors)
        self._transitions = {} if transitions is None else dict(transitions)
        self._shutdown = self.default_shutdown

        self._timeout = timeout

        # Logging to disk.
        self.logging_dir = logging_dir
        self.logger = new_logger(self.logging_dir, str(os.getpid()) + "_main")

    @property
    def phase(self):
        if not self.phases:
            return None

        else:
            return self.phases[0]

    def children_have_obeyed(self, *command_codes):
        if self._child_statuses is None:
            return True

        return all([
            child_status in command_codes
            for child_status in self._child_statuses
        ])

    @property
    def any_children_active(self):
        if self._child_statuses is None:
            return False

        return any([
            child_status == expam.processes.COMMAND_CODES["ACTIVE"]
            for child_status in self._child_statuses
        ])

    @property
    def children_are_dead(self):
        return self.children_have_obeyed(expam.processes.COMMAND_CODES["DEAD"])

    @classmethod
    def from_dict(cls, logging_dir, config):
        def check_int(value):
            if not isinstance(value, int):
                return False

            return True

        if not isinstance(config, dict):
            raise TypeError("Can only be instantiated using dictionary not %s!" % type(config))

        control_center = ControlCenter(logging_dir)

        required_args = [
            "phases",
            "layers",
            "name",
        ]
        optional_args = [
            "timeout",
        ]

        # Set optional arguments first.
        for arg in optional_args:
            if arg in config:

                if arg == "timeout":
                    if type(config[arg]) not in {float, int}:
                        raise ValueError("Parameter must be float or int not %s!" % type(config[arg]))

                    else:
                        control_center._timeout = config["timeout"]

        # Set required arguments.
        for arg in required_args:
            if arg not in config:
                raise ValueError("Required arg %s not in configuration!" % arg)

            else:
                if arg == "phases":
                    control_center.phases = list(config["phases"])

                elif arg == "layers":
                    layers_required_args = {
                        "class",
                        "parent",
                        "child",
                        "n"
                    }
                    current_layer = config["layers"]

                    if not isinstance(current_layer, dict) and "parent" in current_layer:
                        raise ValueError("Incorrect layer configuration!")

                    # Ensure that top layer is directly connected to the control center.
                    if current_layer["parent"] is not ControlCenter:
                        raise TypeError("Top level worker must be connected to ControlCenter not %s!"
                                        % current_layer["parent"])

                    parents = [control_center]

                    at_bottom = False
                    while not at_bottom:

                        if not isinstance(current_layer, dict):
                            raise ValueError("Layer configuration must be dictionary not %s!" % type(current_layer))

                        # Validate arguments.
                        for layer_arg in layers_required_args:
                            if layer_arg not in current_layer:
                                raise ValueError("Missing %s argument required for layer creation!" % layer_arg)

                            elif layer_arg == "n":
                                if not check_int(current_layer["n"]):
                                    raise TypeError("Layer parameter n must be int not %s!" % type(current_layer["n"]))

                            elif layer_arg == "class":
                                if not issubclass(current_layer["class"], expam.processes.JobWorker):
                                    raise TypeError("Process %s must be subclass of JobWorker!"
                                                    % str(current_layer["class"].__name__))

                            elif layer_arg == "child":
                                if current_layer["child"] is None:
                                    at_bottom = True

                        parents = control_center.add_layer(current_layer, parents)
                        current_layer = current_layer["child"]

                elif arg == "name":
                    control_center.group_name = str(config["name"])

        return control_center

    def add_layer(self, current_layer, parents):

        n_children = current_layer["n"]
        n_parents = len(parents)

        def make_shared_array():
            nonlocal n_children

            arr = Array(c_int8, n_children, lock=False)  # Only edited by children.
            for i in range(n_children):
                arr[i] = expam.processes.COMMAND_CODES["ACTIVE"]

            return arr

        def initialise_workers(statuses, timeout):
            nonlocal current_layer, n_children

            # Class of worker.
            worker = current_layer["class"]

            # Custom to initialise the worker.
            if "class_args" in current_layer:

                init_args = [
                    {
                        arg: value
                        if not (isinstance(value, tuple) and len(value) == n_children) else value[cid]
                        for arg, value in current_layer["class_args"].items()
                    }
                    for cid in range(n_children)
                ]

            else:
                init_args = [{} for _ in range(n_children)]

            # Insert default arguments.
            for cid in range(n_children):
                default_args = {
                    "status_index": cid,
                    "statuses": statuses,
                    "child_statuses": None,
                    "timeout": timeout,
                    "phases": self.phases,
                }

                init_args[cid] = {
                    **init_args[cid],
                    **default_args
                }

            return [
                worker(**init_args[cid])
                for cid in range(n_children)
            ]

        def make_pipes():
            nonlocal n_parents, n_children

            parent_cons = [[] for _ in range(n_parents)]
            child_cons = [[] for _ in range(n_children)]

            for i in range(n_parents):
                for j in range(n_children):
                    parent_con, child_con = Pipe()

                    parent_cons[i].append(parent_con)
                    child_cons[j].append(child_con)

            return parent_cons, child_cons

        def make_queues():
            nonlocal n_children

            return [
                Queue()
                for _ in range(n_children)
            ]

        current_statuses = make_shared_array()
        parent_cons, child_cons = make_pipes()
        work_queues = make_queues()

        for i, parent in enumerate(parents):
            parent._child_statuses = current_statuses
            parent._children = parent_cons[i]
            parent.child_queues = work_queues

        workers = initialise_workers(
            current_statuses,
            self._timeout
        )

        for i, worker in enumerate(workers):
            # Set communication routes.
            worker._fellow_statuses = current_statuses
            worker._parents = child_cons[i]
            worker.work_queue = work_queues[i]

        self.workers.append(workers)

        return workers

    def set_methods(self, methods_dict):

        if not isinstance(methods_dict, dict):
            raise TypeError("Config must be dictionary not %s!" % str(type(methods_dict)))

        for key, value in methods_dict.items():
            if value is None:
                continue

            if key == "shutdown":

                if not callable(value):
                    raise TypeError("Shutdown object must be callable!")

                self._shutdown = value

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
                raise ValueError("Unknown key in class config %s!" % key)

    def set_method(self, method_name, phase, method):
        methods = {
            "processor": self._processors,
            "transition": self._transitions,
        }

        if method_name not in methods:
            raise ValueError("Invalid method name %s!" % method_name)

        if method is not None:
            methods[method_name][phase] = method

    def default_shutdown(self):
        self.log("Control center shutting down.")

    def run(self):
        COMMAND_CODES = expam.processes.COMMAND_CODES

        self.start()

        try:

            # Phase loop.
            while self.phases:
                phase = self.phases[0]

                try:
                    self.log("Main filling queue at phase %s..." % self.phase)

                    # Processing loop.
                    while True:
                        # Check children for errors and any custom processing protocols.
                        self.process()

                        for child_queue in self.child_queues:
                            next_job = self.phase_queues[phase].get(block=True, timeout=self._timeout)

                            child_queue.put(next_job)

                except queue.Empty:
                    self.log("Exiting loop from MAIN...")

                except KeyError:

                    # Check there is work for this phase.
                    if phase not in self.phase_queues:
                        pass

                    else:
                        raise

                # Processing finished.
                self.log("Telling children to transition...")
                # Finalise data and transition to next phase.
                self.tell_children(
                    (
                        expam.processes.COMMAND_CODES["SLEEP"],
                        None
                    )
                )

                # Wait for response.
                while not self.children_have_obeyed(COMMAND_CODES["SLEEP"], COMMAND_CODES["DEAD"]):
                    self.process()

                # Transition now that children are transitioning.
                self.transition()  # Progresses current phase.

            # Shutdown phase.
            # Wait for processes to finish.

            while not self.children_have_obeyed(expam.processes.COMMAND_CODES["DEAD"]):
                self.check_children()

                self.timeout()

        except Exception as e:

            # Kill children before raising error.
            # Propagate error downward.
            self.finalise_error(
                -1,
                expam.processes.COMMAND_CODES["ERROR"],
                e
            )

            raise e

        finally:

            self.shutdown()

    def start(self):
        for layer in self.workers:
            for worker in layer:
                worker.start()

    def shutdown(self):
        # Default shutdown procedure.
        self._shutdown()

        # Gather processes.
        self.join()

    def process(self):
        self.check_children()

        if self.phase in self._processors:
            self._processors[self.phase]()

    def transition(self):
        self.log("Main transitioning...")

        if self.phase in self._transitions:
            self._transitions[self.phase]()

        self.phases.pop(0)
        self.log("Main now at phase %s." % ("end" if not self.phases else self.phase))

        # Wake up children and begin next phase.
        self.wake_children()
        self.log("Children are now awake, beginning next phase...")

    def wake_children(self):
        for cid in range(len(self._child_statuses)):
            self._child_statuses[cid] = expam.processes.COMMAND_CODES["ACTIVE"]

    def join(self):
        # Check that my queues are empty.
        all_queues = [*self.child_queues, *[phase_queue for phase_queue in self.phase_queues.values()]]

        for child_queue in all_queues:
            while True:
                try:
                    _ = child_queue.get(block=True, timeout=self._timeout)

                except queue.Empty:
                    break

        # Clean up workers.
        for layer in self.workers:
            for worker in layer:
                worker.join()

    def kill(self):
        for cid in range(len(self._children)):
            self.tell_child(
                cid=cid,
                message=(
                    expam.processes.COMMAND_CODES["DEAD"],
                    None
                )
            )

    def tell_children(self, message):
        for cid in range(len(self._children)):
            self.tell_child(cid, message)

    def tell_child(self, cid, message):
        self.child_queues[cid].put(message)

    def check_children(self):
        for cid in range(len(self._children)):
            self.check_child(cid)

    def check_child(self, cid):
        if self._children[cid].poll(self._timeout):
            err_code, err = self._children[cid].recv()
            self._children[cid].send(True)  # Read receipt.

            self.log("Error received at main!")

            # Propagate error downward.
            self.finalise_error(cid, err_code, err)

            # Processes should have died if an error has been raised.
            self.join()

            print("\n" + "-" * 30 + "\n  Exception raised to MAIN:\n" + "-" * 30)
            raise err

    def finalise_error(self, eid, code, err):
        self.log("MAIN propagating error...")
        self.error = True

        for cid in range(len(self._children)):
            if cid != eid:
                self._children[cid].send((code, err))

        # Wait for children to die.
        while not self.children_are_dead:
            self.timeout()

    def load_job(self, phase, job):
        if phase not in self.phases:
            raise ValueError("Undefined job phase %s!" % phase)

        # Initialise new phase.
        if phase not in self.phase_queues:
            self.phase_queues[phase] = Queue()

        self.phase_queues[phase].put(
            (
                expam.processes.COMMAND_CODES["ACTIVE"],
                job
            )
        )

    def timeout(self):
        time.sleep(self._timeout)

    def log(self, msg):
        self.logger.info(msg)


class ExpamProcesses(ControlCenter):

    def __init__(self, logging_dir, cons, out_dir, *args, **kwargs):
        super(ExpamProcesses, self).__init__(logging_dir, *args, **kwargs)

        # Define non-hierarchical connection to children.
        self.cons = cons
        self.out_dir = out_dir

        # Current UnionWorker array sizes.
        self.arr_sizes = {}

    @classmethod
    def from_method_dict(cls, connections, out_dir, logging_dir, config):

        # Configure base class from configuration dictionary.
        base_center = super(ExpamProcesses, cls).from_dict(logging_dir, config)

        base_arguments = {
            "group_name": "group_name",
            "workers": "workers",
            "child_statuses": "_child_statuses",
            "phases": "phases",
            "phase_queues": "phase_queues",
            "child_queues": "child_queues",
            "children": "_children",
            "processors": "_processors",
            "transitions": "_transitions",
            "timeout": "_timeout",
        }
        base_attributes = {
            attr: getattr(base_center, attr_reference)
            for attr, attr_reference in base_arguments.items()
        }

        # Create child class from this instance.
        control_center = ExpamProcesses(
            logging_dir=logging_dir,
            cons=connections,
            out_dir=out_dir,
            **base_attributes
        )

        # Initialise special methods.
        control_center.set_methods(
            {  # Custom methods for non-hierarchical communication between layers.
                "import": {
                    "processor": control_center.union_checker,
                    "transition": control_center.union_checker,
                },
                "map": {
                    "processor": None,
                    "transition": None,
                },
                "shutdown": control_center.union_saver,
            }
        )

        return control_center

    def union_checker(self):

        for i, con in enumerate(self.cons):
            if con.poll(self._timeout):

                cmd, data = con.recv()

                if cmd == "shm_extend":
                    shm_name, size, create = data

                    if shm_name is None:
                        # Create an allocation.
                        shm_all = shared_memory.SharedMemory(size=size, create=create)

                        shm_name = shm_all.name
                        shm_all.close()

                    else:
                        # Resize a current allocation.
                        shm_all = expam.stores.ResizableSharedMemory(shm_name)

                        shm_all.resize(size)
                        shm_all.close()

                    con.send(shm_name)

                elif cmd == "keys_size":
                    keys_shape = data

                    # Update known size.
                    self.arr_sizes[i] = keys_shape

                    # Set map size.
                    map_shape = keys_shape[:1]

                    # Create a store for its values.
                    shm_values_size = np.product(map_shape) * MAP_DATA_TYPE_SIZE
                    shm_all = shared_memory.SharedMemory(size=shm_values_size, create=True)

                    shm_name = shm_all.name
                    shm_all.close()

                    con.send(shm_name)

                else:
                    raise ValueError("Unrecognised command %s!" % cmd)

    def union_saver(self):

        if not self.error:

            # Calculate combined size.
            final_shape = tuple(
                fn([shape[i] for shape in self.arr_sizes.values()])
                for i, fn in enumerate((sum, max))
            )
            n_cols = final_shape[1]

            # Initialise database on disk.
            db_url = os.path.join(self.out_dir, DB_PATH)
            disk_array = expam.stores.TablesDb(
                db_url=db_url,
                shape=final_shape,
                keys_type=KEYS_DATA_TYPE,
                values_type=MAP_DATA_TYPE,
                create=True
            )
            disk_array.close()

            def add_tuples(*tuples):
                if tuples:
                    return tuple(
                        sum([tup[i] for tup in tuples])
                        for i in range(len(tuples[0]))
                    )

            # Tell each process to write to disk.
            current_shape = (0,) * n_cols

            for con in self.cons:

                while not con.poll(self._timeout):
                    # Don't cause hang due to un-propagated errors!
                    self.check_children()

                cmd, shape = con.recv()

                if cmd != "save":
                    raise ValueError("Unrecognised command %s!" % cmd)

                # Let process know where it comes in the file.
                con.send((current_shape[0], db_url))

                if shape is not None:
                    current_shape = add_tuples(current_shape, shape)

            # Save database parameters to config file.
            save_db_params(
                out_dir=self.out_dir,
                db_params={
                    "keys_shape": current_shape,
                    "keys_data_type": KEYS_DTYPE_STR,
                    "values_shape": current_shape[:1],
                    "values_data_type": MAP_DTYPE_STR,
                }
            )


def main(out_dir, genome_paths, phylogeny_path, k, n=None, n_extract=None,
         n_union=None, pile_size=None, file_name_mode=1, plot=False):
    logs_dir = os.path.join(out_dir, 'logs')

    # Initial time point.
    t0 = datetime.datetime.now()

    # Configure number of processes.
    if n_union is not None and n_extract is not None:
        n_union, n_extract = int(n_union), int(n_extract)

    elif n is not None:
        n_union, n_extract = distribute_processes(n)

    else:
        raise ValueError("Number of processes is can't be defined!")

    # Processing parameter details.
    num_cols = calculate_n_chunks(k)
    ranges = make_ranges(n_union, k)
    genome_paths, maxsize = sort_by_size(genome_paths)

    # Create shm_allocations for in-place kmer processing.
    shm_allocations = prepare_kmer_allocations(maxsize, num_cols, n_extract)

    # Import phylogeny.
    genome_ids, index = import_phylogeny(phylogeny_path)
    # Create LCA matrix.
    format_name = expam.sequences.format_name
    node_to_index = {node.name: i for i, node in enumerate(index.pool) if i > 0}
    lca_matrix = make_lca_matrix(index, node_to_index)

    # Prepare pipes between main and union workers.
    main_cons, union_cons = make_pipes(n_union)
    # Lock system to prevent concurrent access to HDF5 file.
    lock_value = Value(c_uint64, lock=True)

    # Multiprocessing configuration.
    mp_config = {
        "name": "expam",  # Job system title.
        "phases": [  # Two-pass algorithm.
            "import",  # Import sequences and take union.
            "map"  # Map kmer to LCA.
        ],
        "layers": {  # Top layer of child processes - extract sequence & kmers.
            "class": expam.processes.ExtractWorker,
            "class_args": {  # Positional arguments to init class.
                "k": k,
                "shm_kmer_params": shm_allocations,
                "kmer_ranges": ranges,
                "node_map": node_to_index,
                "out_dir": out_dir,
                "logging_dir": logs_dir,
            },
            "parent": ControlCenter,  # Receive work from main process.
            "child": {
                "class": expam.processes.UnionWorker,  # Take union of sets of kmers.
                "class_args": {  # Each either length n_union or 1.
                    "main_con": union_cons,
                    "lca_matrix": lca_matrix,
                    "pile_size": pile_size,
                    "logging_dir": logs_dir,
                    "lock_value": lock_value,
                },
                "parent": expam.processes.ExtractWorker,  # Receive work from ExtractWorkers.
                "child": None,  # Bottom of process hierarchy.
                "n": n_union,  # Number of processes.
            },
            "n": n_extract,
        },
        "timeout": TIMEOUT,  # Time to block for when checking queues.
    }

    process_network = ExpamProcesses.from_method_dict(main_cons, out_dir, logs_dir, mp_config)

    # Load files for importing stage.
    for file_dir in genome_paths:
        process_network.load_job(
            phase="import",
            job=file_dir
        )

    process_network.run()

    # Plot timing data.
    if plot:
        plot_data(logs_dir, t0)

    # Save coordinate names.
    node_names = [node.name for node in index.pool]
    with open(os.path.join(out_dir, COORDS_PATH), "w") as f:
        json.dump(node_names, f)

    # Save LCA matrix.
    np.save(os.path.join(out_dir, LCA_MATRIX_PATH), lca_matrix)

    # Concatenate accession ID files.
    concatenate_accession_ids(os.path.join(out_dir, "phylogeny"))


def save_db_params(out_dir, db_params):
    conf_dir = os.path.join(out_dir, 'expam.conf')

    conf = JSONConfig(url=conf_dir)
    conf.set(**db_params)
    conf.save()


def distribute_processes(n):
    n_union = math.floor(UNION_RATIO * n)
    n_extract = n - n_union

    return n_union, n_extract


def make_ranges(n, k):
    # Calculate hash offset for the k value.
    maxvalue = maximum_hash(k)

    # Interpolate between offset and maxvalue.
    dx = maxvalue // n
    return {
        i: (i * dx, (i + 1) * dx) if i < n - 1
        else (i * dx, maxvalue)
        for i in range(n)
    }


def maximum_hash(k):
    if k < CHUNK_SIZE:
        return 4 ** k - 1
    else:
        return 4 ** CHUNK_SIZE - 1


def sort_by_size(dirs):
    _file_suffixes = expam.sequences.COMP_PARSE
    _suffix_check = expam.sequences.check_suffix

    def _gzip_size(file_dir):
        file_stats = subprocess.run(["gzip", "-l", file_dir], capture_output=True)
        file_meta = [entry for entry in file_stats.stdout.split(b' ') if entry]
        return file_meta[5]

    def _fna_size(file_dir):
        return os.stat(file_dir).st_size

    def _get_file_size(file_dir):
        for suffix in _file_suffixes:
            if _suffix_check(file_dir, suffix):
                return _gzip_size(file_dir)

        return _fna_size(file_dir)

    max_filename = max([len(f) for f in dirs])
    file_dt = np.dtype([
        ("filename", "S%d" % max_filename),
        ("filesize", np.int64)
    ])

    # Turn into structured array.
    file_array = np.array([
        (file_dir, _get_file_size(file_dir))
        for file_dir in dirs
    ], dtype=file_dt)

    # Sort by file size.
    file_array.sort(order="filesize")

    max_size = file_array["filesize"][-1]
    ordered_files = [
        filename.decode("utf8")
        for filename in file_array["filename"].tolist()
    ]

    # Assuming each char in the file is a base.
    return ordered_files, max_size


def prepare_kmer_allocations(rows, cols, n_processes):
    allocation_params = ()
    allocation_size = rows * cols * KEYS_DATA_TYPE_SIZE
    allocation_shape = (rows, cols)

    for _ in range(n_processes):
        # Create and 0 allocation.
        next_shm_allocation = shared_memory.SharedMemory(
            create=True,
            size=allocation_size
        )

        next_arr = np.ndarray(
            shape=allocation_shape,
            dtype=KEYS_DATA_TYPE,
            buffer=next_shm_allocation.buf
        )
        next_arr[:] = 0

        # Pass on allocation details to children.
        allocation_params += ((
                                  next_shm_allocation.name,
                                  allocation_shape
                              ),)

        next_shm_allocation.close()

    return allocation_params


def calculate_n_chunks(k):
    return k // CHUNK_SIZE + (k % CHUNK_SIZE != 0)


def import_phylogeny(phylogeny_path):
    print("Importing phylogeny...")

    # Load the phylogeny into an Index.
    genome_names, index = expam.tree.Index.load_newick(phylogeny_path)
    return genome_names, index


def yield_coordinates(index):
    n = len(index.pool)

    for i in range(1, n):
        for j in range(1, i):
            yield i, j, index.pool[i].coordinate, index.pool[j].coordinate


def compute_lca(i, j, coord_one, coord_two):
    return i, j, expam.tree.propose_lca(coord_one, coord_two)


def make_lca_matrix(index, node_to_index):
    print("Creating LCA matrix...")

    def get_children(fixed_node, flexible_index_list):
        nonlocal index, matrix, node_to_index

        i = 0
        while i < len(flexible_index_list):
            for child in index.pool[flexible_index_list[i]].children:
                m = min(fixed_node, node_to_index[child])
                n = max(fixed_node, node_to_index[child])

                if matrix[n, m] == 0 and n != m:
                    flexible_index_list.append(node_to_index[child])

            i += 1

    n = len(index.pool)
    matrix = np.zeros((n, n), dtype=np.uint16)

    for i in range(1, n):
        subject = index.pool[i].coordinate
        subject_name = index.pool[i].name

        for j in range(1, i):
            if matrix[i, j] > 0:
                continue

            target = index.pool[j].coordinate
            target_name = index.pool[j].name

            lca_coordinate = expam.tree.propose_lca(subject, target)
            lca = index.coord(coordinate=lca_coordinate).name

            left_group = [i]
            right_group = [j]

            if lca not in [subject_name, target_name]:
                # Subject and target are on separate lineages.
                get_children(j, left_group)
                get_children(i, right_group)

            elif lca == subject_name:
                # Target is contained in the clade rooted at subject.
                get_children(i, right_group)

            elif lca == target_name:
                # Subject is contained in the clade rooted at target.
                get_children(j, left_group)

            for p in left_group:
                for q in right_group:
                    m = min(p, q)
                    n = max(p, q)

                    if n != m and matrix[n, m] == 0:
                        matrix[n, m] = node_to_index[lca]

    return matrix


def make_pipes(n):
    parent_cons, child_cons = [], []

    i = 0
    while i < n:
        parent_con, child_con = Pipe()

        parent_cons.append(parent_con)
        child_cons.append(child_con)

        i += 1

    return tuple(parent_cons), tuple(child_cons)


def plot_data(logs_dir: str, t0: datetime.datetime) -> None:
    def _get_time(string):
        nonlocal t0

        date, current_time, am_pm = string.split(" ")[:3]
        mon, day, yr = (int(v) for v in date.split("/"))
        hr, mn, sec = (int(v) for v in current_time.split(":"))

        if "PM" in am_pm and hr < 12:
            hr += 12
        elif "AM" in am_pm and hr == 12:
            hr = 0

        time_point = datetime.datetime(yr, mon, day, hr, mn, sec)

        return (time_point - t0).total_seconds()

    def _get_data(string):
        return float(
            (string.split(" took ")[1])[:-3]
        )

    job_type_map = {  # For labels.
        "request_extension": "remap",
        "disjoint": "disjoint",
        "collapse_pile": "fill",
        "import_sequence": "kmers",
        "send_kmers": "send",
    }
    job_colours = {  # For plotting.
        "request_extension": "red",
        "disjoint": "blue",
        "collapse_pile": "green",
        "import_sequence": "orange",
        "send_kmers": "black",
    }

    process_type_jobs = {
        "_extract": [
            "import_sequence",
            "send_kmers",
        ],
        "_union": [
            "request_extension",
            "disjoint",
            "collapse_pile"
        ],
    }

    fig, ax = plt.subplots(figsize=(20, 15))
    ax.grid()

    log_urls = [
        os.path.join(logs_dir, log_name)
        for log_name in expam.sequences.ls(logs_dir, ext=".log")
        if "_main.log" not in log_name
    ]

    # Plot data one child process at a time.
    for log_url in log_urls:
        process_type = "_extract" if "_extract.log" in log_url else "_union"

        # Data for plotting.
        data = {job_name: ([], []) for job_name in process_type_jobs[process_type]}

        with open(log_url, "r") as f:
            log_data = f.readlines()

        # Read data line by line and append data.
        for i, line in enumerate(log_data):
            for job_name in process_type_jobs[process_type]:

                if job_name in line:
                    data[job_name][0].append(_get_time(line) / 3600)
                    data[job_name][1].append(_get_data(line) / 60)

        # Plot data on same graph.
        for job_name in process_type_jobs[process_type]:
            x, y = np.array(data[job_name][0]), np.array(data[job_name][1])
            ax.plot(
                x,
                y,
                label=job_type_map[job_name],
                marker="x",
                color=job_colours[job_name],
            )

            # Don't allow duplicate labels.
            job_type_map[job_name] = ""

    ax.set(
        xlabel="Time since start [hrs]",
        ylabel="Time taken for task [mins]",
        title="Time taken to do tasks"
    )

    ax.legend(loc='best')
    fig.savefig(os.path.join(logs_dir, "timing.png"))


def timeit(func):
    def inner(*args, **kwargs):
        with Timer() as timer:
            returns = func(*args, **kwargs)

        print(f'@timeit {func.__name__} of {func.__module__}: ', timer)
        return returns

    return inner


def log_timer(func):
    # Time function and log.
    def inner(*args, **kwargs):

        loggers = logging.root.manager.loggerDict
        my_logger = None  # Null case.

        # Appropriate logger will be identified by process id.
        pid = os.getpid()

        for name in loggers:
            if str(pid) in name:
                my_logger = loggers[name]

                break

        if my_logger is None:
            raise Exception("Logger not found!")

        with Timer() as timer:
            returns = func(*args, **kwargs)

        my_logger.info(f'@timeit {func.__name__} of {func.__module__} took {str(timer)}.')
        return returns

    return inner


def new_logger(logs_dir, logger_name):
    """
    Create a fresh Logger instance that doesn't propagate
    to the root logger. Used when a new process is created.

    :param logs_dir:
    :param logger_name:
    :return:
    """
    log_file_url = os.path.join(logs_dir, logger_name + ".log")

    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False  # Don't propagate to root logger.

    if not logger.hasHandlers():  # Don't duplicate logger file handlers.
        # Let each process log to a unique file.
        fh = logging.FileHandler(log_file_url)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(LOG_FORMAT)

        logger.addHandler(fh)

    return logger


def concatenate_accession_ids(f_url):
    def get_data(file_url):
        with open(file_url, "r") as f:
            data = f.read()

        # Delete temporary data file.
        os.remove(file_url)

        return data

    def get_files(folder_url):
        return [
            file_url
            for file_url in expam.sequences.ls(folder_url)
            if "accession_ids_" in file_url
        ]

    id_data = [
        get_data(os.path.join(f_url, file_name))
        for file_name in get_files(f_url)
    ]

    with open(os.path.join(f_url, "accession_ids.csv"), "w") as f:
        f.write("\n".join(id_data))
