from multiprocessing import Process
import queue
import time
import traceback

import psutil

from expam.database import TIMEOUT
from expam.logger import new_logger
from expam.process import COMMAND_CODES


class JobWorker(Process):
    def __init__(
        self, logging_dir, status_index, statuses,
        child_statuses, timeout=TIMEOUT, job_name="Worker",
        phases=None
    ):
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

