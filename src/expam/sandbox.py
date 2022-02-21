#!/usr/bin/env python3

import logging
import multiprocessing as mp
import psutil
import re
import resource
import subprocess
import time


def limit(mem_limit):
    """
    Defines a wrapper for functions that first sets a resource limit before
    executing the function.

    :param mem_limit:
    :return:
    """
    def inner(fn):
        """
        Wrapper that sets some memory limit then executes the input function.

        :param fn:
        :return:
        """
        def limited_fn(*args, **kwargs):
            soft, hard = resource.getrlimit(resource.RLIMIT_RSS)
            resource.setrlimit(resource.RLIMIT_RSS, (mem_limit, hard))

            fn(*args, **kwargs)

        return limited_fn
    return inner


def sandbox(logging_name, mem_limit, interval):
    """
    Create a decorator that runs function in a new process, limiting the resources of that
    process and observing it at regular intervals.

    :param mem_limit:
    :param interval:
    :return:
    """
    def inner(fn):
        """
        Run function in a new process, modifying the function so that it first sets a resource
        limit for the environment it is in.

        :param fn:
        :return:
        """
        def job(*args, **kwargs):
            # Configure logging output.
            logger = logging.getLogger(logging_name)
            logger.info("\ttotal\tused\tfree\tshared\tbuff/cache\tavailable")

            # Start new child process.
            ctx = mp.get_context('fork')

            new_target = limit(mem_limit)(fn)
            p = ctx.Process(target=new_target, args=args, kwargs=kwargs)

            p.start()

            while p.is_alive():
                # Log memory usage.
                mem_info = subprocess.run(['free', '-h'], capture_output=True).stdout.decode()
                logger.info("\t".join(re.findall(r"[^\s]+", mem_info)[6:13]))

                time.sleep(interval)

            p.join()

        return job
    return inner


def parse(cmd, mem_limit=None, interval=1, out=None):
    i = 0

    if len(cmd) > 0 and 'sandbox' in cmd[i]:
        i += 1

    while i < len(cmd):
        if cmd[i] in {'--memory', '-m'}:
            mem_limit = int(cmd[i + 1])
            i += 1

        elif cmd[i] in {'--interval', '-t'}:
            interval = float(cmd[i + 1])
            i += 1

        elif cmd[i] in {'--out', '-o'}:
            out = cmd[i + 1]
            i += 1

        elif cmd[i] in {'--x', '-x'}:
            total_memory = psutil.virtual_memory().total
            mem_limit = int(float(cmd[i + 1]) * total_memory)
            i += 1

        else:
            break

        i += 1

    cmd = cmd[i:]
    return mem_limit, interval, out, cmd


def sandbox_command(data):
    mem_limit, interval, out_file, cmd = parse(data)

    if mem_limit is None:
        soft, hard = resource.getrlimit(resource.RLIMIT_RSS)
        mem_limit = soft

    logger = logging.getLogger('sandbox_application')
    logger.setLevel(logging.DEBUG)

    # Create the appropriate handler.
    if out_file is None:
        handle = logging.StreamHandler()
    else:
        handle = logging.FileHandler(out_file)

    handle.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s ... %(message)s')
    handle.setFormatter(formatter)

    logger.addHandler(handle)

    @sandbox(logging_name='sandbox_application', mem_limit=mem_limit, interval=interval)
    def run_cmd(command):
        subprocess.run(command, shell=False, check=True)

    if cmd:
        run_cmd(cmd)


def main():
    import sys
    sandbox_command(sys.argv)


if __name__ == '__main__':
    main()
