from argparse import ArgumentParser, Namespace, RawTextHelpFormatter
from collections import namedtuple
import datetime
import os
import shutil
from typing import Set

import matplotlib.pyplot as plt
import numpy as np

from expam import __version__
from expam.utils import die, ls, make_path_absolute, parse_float, parse_int

ExpamOptions = namedtuple(
    'ExpamOptions',
    [
        # Runtime arguments
        'command', 'db_name', 'k', 'n', 's', 'phylogeny', 'alpha',
        # Directory arguments
        'directory', 'out_url', 'truth_dir',
        # Parameter arguments
        'length', 'pile', 'error_rate', 'first_n', 'paired_end',
        # Summary arguments
        'plot', 'cutoff', 'cpm', 'taxonomy',
        # Plot arguments
        'groups', 'phyla', 'keep_zeros', 'ignore_names', 'colour_list', 'rank', 'log_scores', 'itol_mode',
        # Tree arguments
        'use_sourmash', 'use_rapidnj', 'use_quicktree'
    ]
)

def retrieve_arguments() -> ExpamOptions:
    parser = ArgumentParser(description="  expam CLI\n--------------\n", formatter_class=RawTextHelpFormatter)

    parser.add_argument('--version', action='version', version='%(prog)s ' + ".".join(str(v) for v in __version__))
    parser.add_argument("command", default=None,
                        help='\nCommand to execute. Valid commands include:\n'
                             '-------------------------------------------\n'
                             'create:-\tInitialise database.\n'
                             'build:-\t\tStart building database.\n'
                             'print:-\t\tPrint current database parameters.\n'
                             'run:-\t\tRun reads against database.\n'
                             'add:-\t\tAdd sequence to the database.\n'
                             'remove:-\tRemove sequence from database (only impacts future db builds).\n'
                             'set:-\t\tSet database build parameters.\n'
                             'to_taxonomy:-\t\tConvert results to taxonomic setting.\n'
                             'phylotree:-\t\tDraw results on phylotree.\n'
                             'draw_tree:-\t\tDraw the reference tree.\n'
                             'download_taxonomy:-\t\tDownload taxonomic information for reference seqeunces.\n'
                             'cutoff:-\t\tApply cutoff to some set of already processed classifications. THIS WILL OVERWRITE OLD RESULTS!\n'
                             'mashtree:-\tCreate mashtree from current sequences and add to database.\n'
                             'quickrun:-\tInitialise, set parameters and start building db (assumes\n'
                             '\t\t\tsequences all lie in the same folder).\n'
                             'make_reads:-\tUniformly sample reads of length l from some input sequence.\n'
                             '\t\tThis is for testing purposes only, and is not a replacement\n'
                             '\t\tfor actual read generating software.\n',
                        metavar="[command]")
    parser.add_argument("-db", "--db_name", dest="db_name",
                        help="Name of database.",
                        metavar="[database name]")
    parser.add_argument("-k", "--kmer", dest="k",
                        help="Length of mer used for analysis.",
                        metavar="[k value (int)]")
    parser.add_argument("-n", "--n-processes", dest="n",
                        help="Number of CPUs to use for processing.",
                        metavar="[n (int)]")
    parser.add_argument("-s", "--sketch", dest="s",
                        help="Sketch size for mash.",
                        metavar="[sketch size (int)]")
    parser.add_argument("-p", "--phylogeny", dest="phylogeny",
                        help="URL of Newick file containing phylogeny.",
                        metavar="[phylogeny URL]")
    parser.add_argument("-d", "--directory", dest="directory", action="append",
                        help="File URL, context depending on command supplied.",
                        metavar="[directory]")
    parser.add_argument("-l", "--length", dest="length",
                        help="Length of simulated reads.",
                        metavar="[read length]")
    parser.add_argument("-o", "--out", dest="out_url",
                        help="Where to save classification results.",
                        metavar="[out URL]")
    parser.add_argument("-y", "--pile", dest="pile",
                        help="Number of genomes to pile at a time (or inf).",
                        metavar="[pile size]")
    parser.add_argument("-e", "--error-rate", dest="error_rate",
                        help="Generate error in reads (error ~ reads with errors / reads).",
                        metavar="[error rate]")
    parser.add_argument("-t", "--truth", dest="truth_dir",
                        help="Location of truth dataset.")
    parser.add_argument("--plot", dest="plot", default=False, action="store_true",
                        help="Plot timing data of database build.")
    parser.add_argument("--first", dest="first_n", default=None,
                        help="Add first n genomes in folder.")
    parser.add_argument("--cutoff", dest="cutoff", default=0,
                        help="Ignore organisms with less than `cutoff` reads in results.")
    parser.add_argument("--cpm", dest="cpm", default=100,
                        help="Counts/million cutoff for read-count to be non-negligible.")
    parser.add_argument("--taxonomy", dest="taxonomy", default=False, action="store_true",
                        help="Convert phylogenetic results to taxonomic results.")
    parser.add_argument("--phyla", dest="phyla", default=False, action="store_true",
                        help="Colour phylotree results by phyla.")
    parser.add_argument("--rank", dest="rank", default=None,
                        help="Rank at which to sort results.")
    parser.add_argument("--keep-zeros", dest="keep_zeros", default=False, action="store_true",
                        help="Keep nodes of output where no reads have been assigned.")
    parser.add_argument("--ignore-names", dest="ignore_names", default=False, action="store_true")
    parser.add_argument("--group", dest="groups", action="append", nargs="+",
                        help="Space-separated list of sample files to be treated as a single group in phylotree.")
    parser.add_argument("--colour-list", dest="colour_list", nargs="+",
                        help="List of colours to use when plotting groups in phylotree.")
    parser.add_argument("--sourmash", dest="use_sourmash", default=False, action="store_true",
                        help="Use sourmash for distance estimation.")
    parser.add_argument("--rapidnj", dest="use_rapidnj", default=True, action="store_true",
                        help="Use RapidNJ for Neighbour-Joining algorithm.")
    parser.add_argument("--quicktree", dest="use_quicktree", default=False, action="store_true",
                        help="Use QuickTree for Neighbour-Joining algorithm.")
    parser.add_argument("--paired", dest="paired_end", default=False, action="store_true",
                        help="Treat reads as paired-end.")
    parser.add_argument("--alpha", dest="alpha", default=0.1,
                        help="Percentage requirement for classification subtrees (see Tutorials 1 & 2).")
    parser.add_argument("--log-scores", dest="log_scores", default=False, action="store_true",
                        help="Log transformation to opacity scores on phylotree (think uneven distributions).")
    parser.add_argument("--itol", dest="itol_mode", default=False, action="store_true",
                        help="Output plotting data in ITOL format.")

    # Parse arguments.
    args: Namespace = parser.parse_args()

    if args.db_name is None:
        args.db_name = db_name_from_environment()
    else:
        args.db_name = make_path_absolute(args.db_name, os.getcwd())

    return ExpamOptions(**{field: getattr(args, field) for field in ExpamOptions._fields})


def db_name_from_environment():
    try:
        return os.environ["EXPAM_DB_DIR"]
    except KeyError:
        return None


def clear_logs(log_path) -> None:
    print("Clearing old log files...")

    try:
        shutil.rmtree(log_path)
    except FileNotFoundError:
        pass

    try:
        os.mkdir(log_path)
    except OSError:
        die("Can't make log path %s!" % log_path)


class CommandGroup:
    commands: Set[str] = {}

    @classmethod
    def take_args(cls, args: ExpamOptions) -> dict:
        return {}

    def run(self, command) -> None:
        try:
            getattr(self, command)()
        except AttributeError:
            raise AttributeError("Command %s not found!" % command)

    @classmethod
    def parse_ints(cls, *params):
        if len(params) == 1:
            return parse_int(params[0])
        else:
            return (parse_int(param) for param in params)

    @classmethod
    def parse_floats(cls, *params):
        if len(params) == 1:
            return parse_float(params[0])
        else:
            return (parse_float(param) for param in params)

    @staticmethod
    def get_user_confirmation(msg):
        if "y/n" not in msg:
            msg += " (y/n)"

        while True:
            ans = input(msg)
            ans = ans.lower().strip()

            if ans == "y":
                return True
            elif ans == "n":
                return False
            else:
                die("Invalid response, please provide y/n.")

    @staticmethod
    def get_time(data) -> datetime.datetime:
        date_data, time_data, am_pm = data.split(" ")[:3]

        mon, day, year = (int(v) for v in date_data.split("/"))
        hr, min, sec = (int(v) for v in time_data.split(":"))

        if "PM" in am_pm and hr < 12:
            hr += 12
        elif "AM" in am_pm and hr == 12:
            hr = 0

        return datetime.datetime(year, mon, day, hr, min, sec)

    def get_t0(self, logs_dir: str):
        def _first(itr):
            try:
                return itr[0]
            except IndexError:
                raise ValueError("Can't find main log file!")

        # Get effective t0 from main log.
        main_log_dir = _first([
            os.path.join(logs_dir, log_name)
            for log_name in ls(logs_dir, ext=".log")
            if "_main" in log_name
        ])

        with open(main_log_dir, "r") as f:
            log_data = f.readline()  # Get time of first log point.

        return main_log_dir, self.get_time(log_data)


class PlotLogs:
    def __init__(self, logs_dir: str, t0: datetime.datetime) -> None:
        self.logs_dir = logs_dir
        self.t0 = t0

    def plot(self):
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
            os.path.join(self.logs_dir, log_name)
            for log_name in ls(self.logs_dir, ext=".log")
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
                        data[job_name][0].append(self._get_time(line) / 3600)
                        data[job_name][1].append(self._get_data(line) / 60)

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
        fig.savefig(os.path.join(self.logs_dir, "timing.png"))

    def _get_time(self, string):
        date, current_time, am_pm = string.split(" ")[:3]
        mon, day, yr = (int(v) for v in date.split("/"))
        hr, mn, sec = (int(v) for v in current_time.split(":"))

        if "PM" in am_pm and hr < 12:
            hr += 12
        elif "AM" in am_pm and hr == 12:
            hr = 0

        time_point = datetime.datetime(yr, mon, day, hr, mn, sec)

        return (time_point - self.t0).total_seconds()

    @staticmethod
    def _get_data(string):
        return float(
            (string.split(" took ")[1])[:-3]
        )
