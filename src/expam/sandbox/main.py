#!/usr/bin/env python3

import logging
import multiprocessing as mp
import sys
import psutil
import re
import resource
import subprocess
import time


def check_free_command():
    err = subprocess.run(['free', '--mega'], capture_output=True, shell=True).stderr

    if err:
        return False
    else:
        return True


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

            # Initial log for starting memory usage.
            mem_info = subprocess.run(['free', '--mega'], capture_output=True).stdout.decode()
            logger.info("Initial...%s" % "\t".join(re.findall(r"[^\s]+", mem_info)[6:13]))

            # Start new child process.
            ctx = mp.get_context('fork')

            new_target = limit(mem_limit)(fn)
            p = ctx.Process(target=new_target, args=args, kwargs=kwargs)

            p.start()

            while p.is_alive():
                # Log memory usage.
                mem_info = subprocess.run(['free', '--mega'], capture_output=True).stdout.decode()
                logger.info("\t".join(re.findall(r"[^\s]+", mem_info)[6:13]))

                time.sleep(interval)

            p.join()

        return job
    return inner


def parse():
    from argparse import ArgumentParser, RawTextHelpFormatter, REMAINDER

    parser = ArgumentParser(description="::: expam_limit ::: \n\n\tcommand for limiting and logging memory usage.", formatter_class=RawTextHelpFormatter)
    parser.add_argument('-m', '--memory', type=int, dest='memory', help="Explicit memory limit in bytes.")
    parser.add_argument('-x', '--x', type=float, dest='x', help="Memory limit as a percentage of total system memory.")
    parser.add_argument('-t', '--interval', type=float, dest='interval', default=1.0, help="Time between memory logs.")
    parser.add_argument('-o', '--out', type=str, dest='out', help="File to log to (default is stdout).")
    parser.add_argument('command', nargs=REMAINDER, help="Command to run in limited regime.")

    args = parser.parse_args()

    # Establish memory limit.
    byte_limit = args.memory
    perc_limit = None if args.x is None else int(args.x * psutil.virtual_memory().total)

    if byte_limit is not None and perc_limit is not None:
        mem_limit = min(byte_limit, perc_limit)

    elif byte_limit is not None:
        mem_limit = byte_limit
    
    elif perc_limit is not None:
        mem_limit = byte_limit
    
    else:
        mem_limit = None

    return mem_limit, args.interval, args.out, args.command


def sandbox_command():
    mem_limit, interval, out_file, cmd = parse()

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
    if not check_free_command():
        print("Must run on Linux system (or have `free` command available).")
        sys.exit(1)
    
    sandbox_command()


if __name__ == '__main__':
    main()
