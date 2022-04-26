import multiprocessing
import platform
from typing import Tuple

from expam.cli.build import BuildCommand
from expam.cli.classify import ClassifyCommand
from expam.cli.main import CommandGroup, ExpamOptions, retrieve_arguments
from expam.cli.tree import TreeCommand
from expam.cli.utils import UtilsCommand
from expam.utils import die


def main():
    # Set process creation method.
    if platform.system() != "Windows":
        multiprocessing.set_start_method("fork")

    args: ExpamOptions = retrieve_arguments()

    handlers: Tuple[CommandGroup] = (BuildCommand, ClassifyCommand, TreeCommand, UtilsCommand)
    for handler in handlers:
        if args.command in handler.commands:
            handler(**handler.take_args(args)).run(args.command)
            break
    else:
        die("Unrecognised command: %s" % args.command)


if __name__ == "__main__":
    main()
        