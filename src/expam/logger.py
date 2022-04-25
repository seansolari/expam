import logging
import os
from timeit import default_timer
from expam.database import LOG_FORMAT

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
