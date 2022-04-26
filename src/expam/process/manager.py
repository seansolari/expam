from ctypes import *
from multiprocessing import Array, Pipe, Queue, shared_memory
import os
import queue
import time

import numpy as np
from expam.database import TIMEOUT, DataTypeConfig, FileLocationConfig, expam_dtypes
from expam.database.config import JSONConfig, load_database_config
from expam.database.db import TablesDb
from expam.database.shm import ResizableSharedMemory
from expam.logger import new_logger
from expam.process import COMMAND_CODES
from expam.process.jobworker import JobWorker


class ControlCenter:
    def __init__(
        self, logging_dir, group_name=None, workers=None, child_statuses=None, phases=None,
        phase_queues=None, child_queues=None, children=None, processors=None, transitions=None,
        timeout=TIMEOUT
):
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
            child_status == COMMAND_CODES["ACTIVE"]
            for child_status in self._child_statuses
        ])

    @property
    def children_are_dead(self):
        return self.children_have_obeyed(COMMAND_CODES["DEAD"])

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
                                if not issubclass(current_layer["class"], JobWorker):
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
                arr[i] = COMMAND_CODES["ACTIVE"]

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
                        COMMAND_CODES["SLEEP"],
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

            while not self.children_have_obeyed(COMMAND_CODES["DEAD"]):
                self.check_children()

                self.timeout()

        except Exception as e:

            # Kill children before raising error.
            # Propagate error downward.
            self.finalise_error(
                -1,
                COMMAND_CODES["ERROR"],
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
            self._child_statuses[cid] = COMMAND_CODES["ACTIVE"]

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
                    COMMAND_CODES["DEAD"],
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
                COMMAND_CODES["ACTIVE"],
                job
            )
        )

    def timeout(self):
        time.sleep(self._timeout)

    def log(self, msg):
        self.logger.info(msg)


class ExpamProcesses(ControlCenter):
    def __init__(self, logging_dir, cons, out_dir, dtype_config: DataTypeConfig = expam_dtypes, *args, **kwargs):
        super(ExpamProcesses, self).__init__(logging_dir, *args, **kwargs)

        # Define non-hierarchical connection to children.
        self.cons = cons
        self.out_dir = out_dir

        self.dtypes: DataTypeConfig = dtype_config
        self.database_config: FileLocationConfig = load_database_config(out_dir)

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
                        shm_all = ResizableSharedMemory(shm_name)

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
                    shm_values_size = np.product(map_shape) * self.dtypes.map_dtype_size
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
            db_url = self.database_config.database_file
            disk_array = TablesDb(
                db_url=db_url,
                shape=final_shape,
                keys_type=self.dtypes.keys_dtype,
                values_type=self.dtypes.map_dtype,
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
            self.save_db_params(
                db_params={
                    "keys_shape": current_shape,
                    "keys_data_type": self.dtypes.keys_dtype_str,
                    "values_shape": current_shape[:1],
                    "values_data_type": self.dtypes.map_dtype_str,
                }
            )

    def save_db_params(self, db_params):
        conf = JSONConfig(url=self.database_config.conf)
        conf.set(**db_params)
        conf.save()

