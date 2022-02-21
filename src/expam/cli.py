import datetime
from functools import partial
import json.decoder
import math
import platform
import re

import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
from argparse import ArgumentParser, RawTextHelpFormatter

from expam.run import die

EXPAM_VERSION = (0, 0, 7)
print("expam v%s" % ".".join(str(i) for i in EXPAM_VERSION))

# Check python version.
if sys.version_info[0] < 3:
    raise Exception("Must be using Python3.8 or later.")
elif sys.version_info[1] < 8:
    raise Exception("Must be using Python3.8 or later.")

# Determine current directory.
CWD = os.getcwd()

# Commands that don't require db_name parameter.
AUTO_COMMANDS = {
    "count",
    "from_jellyfish",
    "make_reads",
    "make_phylogeny",
}


#
#   Extract datetime point from string.
#
def get_time(data):
    date_data, time_data, am_pm = data.split(" ")[:3]

    mon, day, year = (int(v) for v in date_data.split("/"))
    hr, min, sec = (int(v) for v in time_data.split(":"))

    if "PM" in am_pm and hr < 12:
        hr += 12
    elif "AM" in am_pm and hr == 12:
        hr = 0

    return datetime.datetime(year, mon, day, hr, min, sec)


#
# Find the time of the first logging message from a set of log files at logs_dir.
#
def get_t0(logs_dir):
    from expam.sequences import ls

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

    return main_log_dir, get_time(log_data)


#
# Request boolean confirmation from user.
#
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


#
# Check that supplied parameters are integers, and if so
# return integer form.
#
def check_int_params(*params):
    for param in params:
        INVALID_PARAM_MSG = ("Invalid parameter (%s), must be integer!" % str(param))

        if param is not None:

            try:
                # Convert to float.
                param = float(param)

            except ValueError:
                die(INVALID_PARAM_MSG)

            # Convert to int and see if the value changes.
            new_param = int(param)
            if new_param != param:
                die(INVALID_PARAM_MSG)

            param = new_param

        yield param


#
# Try to return float form of arguments.
#
def check_float_params(*params):
    for param in params:
        INVALID_PARAM_MSG = ("Invalid parameter (%s), must be integer!" % str(param))

        if param is not None:
            try:
                param = float(param)

            except ValueError:
                die(INVALID_PARAM_MSG)

        yield param


def check_phylogeny(conf):
    if conf['phylogeny_path'] is None:
        return False
    return True


#
# Load a database configuration file.
#
def load_configuration_file(out_dir, return_config=False):
    from expam.run import JSONConfig

    # Load config file.
    try:
        configParser = JSONConfig(url=os.path.join(out_dir, "expam.conf"))
    except json.decoder.JSONDecodeError:
        die("Invalid configuration file!")
        return

    if return_config:
        return configParser

    genome_paths, phylogeny_path, k, n, pile_size = (
        configParser.get_paths(), configParser["phylogeny_path"],
        int(configParser["k"]), int(configParser["n"]), configParser["pile"])

    # Format phylogeny path.
    phylogeny_path = make_path_absolute(phylogeny_path, out_dir)

    return k, n, phylogeny_path, genome_paths, pile_size


#
# Calculate alpha to allow for at most one SNP. Assume a read length of 150.
#
def calculate_alpha(k):
    return 1e-3 * math.floor(((150 - k) / 150) * 1e3)


def _checker(cmd):
    _has_mash = shutil.which(cmd) is not None
    if not _has_mash:
        die('Could not find local installation of %s!' % cmd)


def _check_mash():
    _checker('mash')


def _check_mashtree():
    _checker('mashtree')


def _check_rapidnj():
    _checker('rapidnj')


def _check_quicktree():
    _checker('quicktree')


def subp(command, cwd=CWD):
    p = subprocess.run(
        command,
        shell=True,
        capture_output=True,
        cwd=cwd
    )

    # Check for fails.
    if p.returncode != 0:
        raise Exception("\n" + p.stderr.decode("utf-8"))

    return p.stdout.decode("utf-8")


def mash_sketch(k, s, p, sequences, out_dir):
    _check_mash()

    cmd_fmt = "mash sketch -k %d -p %d -s %d -o %s %s"
    cmd = cmd_fmt % (
        k,
        p,
        s,
        out_dir,
        ' '.join(sequences)
    )

    subp(cmd)


def sour_sketch(k, s, n, sequences, sig_dir):
    from expam.sourtree import make_signatures

    make_signatures(n, sequences, sig_dir, k, s)


def get_groups(conf, group_input):
    groups = conf.groups() if group_input is None else [group_input]

    # Check group has more than zero sequences.
    for group in groups:
        k, s, sequences = conf.group_get(group)

        if len(sequences) > 0:
            if k is None or s is None:
                raise ValueError("Parameters unspecified for group %s (k=%s, s=%s)!" % (group, str(k), str(s)))

            yield group


def check_sketches(conf, phy_dir, group, use_sourmash=False):
    sketch_dir = os.path.join(phy_dir, 'sketch')
    sketch_name_fmt = "%s.k%d.s%d.%s"

    if not os.path.exists(sketch_dir):
        return False

    for group_name in get_groups(conf, group):
        k, s, sequences = conf.group_get(group_name)
        file_name = sketch_name_fmt % (group_name, k, s, "%s")

        dest = os.path.join(sketch_dir, file_name % ("sour" if use_sourmash else "msh"))
        if not os.path.exists(dest):
            return False

    return True


def do_sketches(conf, phy_dir, group, use_sourmash=False):
    print("Sketching sequences...")

    sketch_dir = os.path.join(phy_dir, 'sketch')
    n_procs = conf['n'] if conf['n'] is not None else 1

    if not os.path.exists(sketch_dir):
        os.mkdir(sketch_dir)

    sketch_name_fmt = "%s.k%d.s%d.%s"

    for group_name in get_groups(conf, group):
        k, s, sequences = conf.group_get(group_name)
        file_name = sketch_name_fmt % (group_name, k, s, "%s")

        if use_sourmash:
            sig_path = os.path.join(sketch_dir, file_name % "sour")

            if not os.path.exists(sig_path):
                os.mkdir(sig_path)

            sour_sketch(k=k, s=s, n=n_procs, sequences=sequences, sig_dir=sig_path)

        else:
            _check_mash()

            # Use mash by default.
            file_path = os.path.join(sketch_dir, file_name % "msh")
            mash_sketch(k=k, s=s, p=n_procs, sequences=sequences, out_dir=file_path)


def check_distances(conf, phy_dir, group):
    sketch_dir = os.path.join(phy_dir, 'sketch')
    dist_dir = os.path.join(phy_dir, 'distance')

    if not os.path.exists(dist_dir):
        return False

    dist_name_fmt = "%s.k%d.s%d.tab"

    for group_name in get_groups(conf, group):
        k, s, sequences = conf.group_get(group_name)
        file_name = dist_name_fmt % (group_name, k, s)

        dest = os.path.join(sketch_dir, file_name)
        if not os.path.exists(dest):
            return False

    return True


def do_distances(conf, phy_dir, group, use_sourmash=False):
    print("Calculating pairwise distances...")

    sketch_dir = os.path.join(phy_dir, 'sketch')
    dist_dir = os.path.join(phy_dir, 'distance')
    n_procs = conf['n'] if conf['n'] is not None else 1

    if not os.path.exists(dist_dir):
        os.mkdir(dist_dir)

    dist_name_fmt = "%s.k%d.s%d.%s"

    for group_name in get_groups(conf, group):
        k, s, sequence_files = conf.group_get(group_name)

        sketch_name = dist_name_fmt % (group_name, k, s, '%s')

        matrix_name = dist_name_fmt % (group_name, k, s, 'tab')
        matrix_file = os.path.join(dist_dir, matrix_name)

        if use_sourmash:
            sketch_dir = os.path.join(sketch_dir, sketch_name % 'sour')
            sour_dist(n=n_procs, sig_dir=sketch_dir, matrix_dir=matrix_file)

        else:
            _check_mash()

            sketch_file = os.path.join(sketch_dir, sketch_name % 'msh')
            mash_dist(n_procs, sketch_file, matrix_file)


def check_trees(conf, phy_dir, group):
    tree_dir = os.path.join(phy_dir, 'tree')

    if not os.path.exists(tree_dir):
        return False

    for group_name in get_groups(conf, group):
        tree_name = '%s.nwk' % group_name
        tree_path = os.path.join(tree_dir, tree_name)

        if not os.path.exists(tree_path):
            return False

    return True


def do_trees(conf, phy_dir, group, use_quicktree=False):
    print("Running neighbour-joining...")

    dist_dir = os.path.join(phy_dir, 'distance')
    tree_dir = os.path.join(phy_dir, 'tree')

    n_procs = conf['n'] if conf['n'] is not None else 1
    dist_fmt = "%s.k%d.s%d.tab"

    if not os.path.exists(tree_dir):
        os.mkdir(tree_dir)

    for group_name in get_groups(conf, group):
        k, s, _ = conf.group_get(group_name)

        matrix_name = dist_fmt % (group_name, k, s)
        matrix_path = os.path.join(dist_dir, matrix_name)

        tree_name = "%s.nwk" % group_name
        tree_path = os.path.join(tree_dir, tree_name)

        if use_quicktree:
            quicktree(matrix_dir=matrix_path, tree_dir=tree_path)

        else:
            rapidnj(matrix_dir=matrix_path, n=n_procs, tree_dir=tree_path)


def mash_dist(p, sketch_dir, matrix_dir):
    _check_mash()

    cmd_fmt = "mash dist -p %d -t %s %s"
    cmd = cmd_fmt % (
        p,
        sketch_dir,
        sketch_dir
    )

    _unfmt_dst = subp(cmd)
    lines = [line.strip() for line in _unfmt_dst.split('\n')]

    names = lines[0].split('\t')
    lines[0] = '%d' % len(names[1:])

    # Format sequences names from within directories.
    for i in range(1, len(lines)):
        line = lines[i].split('\t')
        line[0] = os.path.basename(line[0])
        lines[i] = '\t'.join(line)

    with open(matrix_dir, 'w') as f:
        f.write('\n'.join(lines))


def sour_dist(n, sig_dir, matrix_dir):
    from expam.sourtree import make_distances

    make_distances(sig_dir, n, matrix_dir)


def _format_tree_string(nwk):
    seq_types = [".fna", ".faa"]
    comp_types = [".tar.gz", ".gz"]

    nwk = nwk.replace("'", "")

    for seq_type in seq_types:
        for comp_type in comp_types:
            ext = seq_type + comp_type + ":"
            nwk = nwk.replace(ext, ":")

    return nwk


def rapidnj(matrix_dir, n, tree_dir):
    _check_rapidnj()

    cmd_fmt = "rapidnj %s -i pd -o t -c %d"
    cmd = cmd_fmt % (matrix_dir, n)

    unformatted_tree = subp(cmd)
    tree = _format_tree_string(unformatted_tree)

    with open(tree_dir, 'w') as f:
        f.write(tree)

    print("RapidNJ wrote tree to %s." % tree_dir)


def quicktree(matrix_dir, tree_dir):
    _check_quicktree()

    cmd_fmt = "quicktree -in m -out t %s"
    cmd = cmd_fmt % matrix_dir

    unformatted_tree = subp(cmd)
    tree = _format_tree_string(unformatted_tree)

    with open(tree_dir, 'w') as f:
        f.write(tree)

    print("QuickTree wrote tree to %s." % tree_dir)


def do_mashtree(conf, phy_dir, group):
    print("Creating mashtree...")
    _check_mashtree()

    tree_dir = os.path.join(phy_dir, 'tree')
    temp_dir = os.path.join(phy_dir, 'tmp')

    n_procs = conf['n'] if conf['n'] is not None else 1

    if not os.path.exists(tree_dir):
        os.mkdir(tree_dir)

    for group_name in get_groups(conf, group):
        k, s, sequences = conf.group_get(group_name)
        tree_path = os.path.join(tree_dir, "%s.nwk" % group_name)

        mashtree(k, n_procs, s, sequences, tree_path, temp_dir)


def mashtree(k, n, s, sequences, tree_dir, temp_dir):
    _names_file = os.path.join(temp_dir, 'sequence_names.txt')

    def delete_temp():
        nonlocal temp_dir

        shutil.rmtree(temp_dir)
        print("Deleted temporary directory %s." % temp_dir)

    # Make temporary directory.
    try:
        print("Making temporary directory %s..." % temp_dir)
        os.mkdir(temp_dir)

    except OSError:
        try:
            delete_temp()
        except FileNotFoundError:
            pass

        die("Could not create temporary directory!")

    # Write file of files.
    try:
        with open(_names_file, "w") as f:
            f.write("\n".join(sequences))

        # Make mashtree and insert into configuration file.
        print("Making mashtree...")

    except OSError:
        die("Could not write genome paths to text file!")

    cmd = "mashtree --numcpus %d --kmerlength %d --sketch-size %d --file-of-files %s" \
          % (n, k, s, _names_file)

    # Build the tree.
    try:
        return_val = subp(cmd, cwd=temp_dir)

        with open(tree_dir, "w") as f:
            f.write(return_val)

        print("Phylogeny created in %s." % tree_dir)

    except:
        try:
            delete_temp()
        except FileNotFoundError:
            pass

        raise

    # Delete the temporary directory.
    try:
        delete_temp()

    except OSError:
        die("Could not delete temporary directory %s!" % temp_dir)


def finalise_tree(conf, phy_dir, db_name):
    print("Finalising tree...")

    from expam.sequences import ls

    tree_dir = os.path.join(phy_dir, 'tree')
    tree_name = "%s.nwk" % db_name
    tree_path = os.path.join(tree_dir, tree_name)

    tree_files = [os.path.join(tree_dir, name) for name in ls(tree_dir, ext='.nwk') if name != tree_name]

    def get_tree_data(file_url):
        with open(file_url, 'r') as f:
            tree_data = f.read().strip().rstrip(";")

        return tree_data

    if len(tree_files) == 1:
        nwk_data = get_tree_data(tree_files[0]) + ";"

        with open(tree_path, 'w') as f:
            f.write(nwk_data)

    else:
        try:
            template = get_tree_data(tree_path) + ";"

        except FileNotFoundError:
            raise Exception("No template found! Please write a template to %s!" % tree_path)

        template_groups = re.findall(r"{{(\S+?)}}", template)

        tree_data = {
            group_name: get_tree_data(tree_file)
            for group_name in conf.groups()
            for tree_file in tree_files
            if group_name in tree_file
        }

        for template_group in template_groups:
            if template_group not in tree_data:
                die("Template group %s has no tree!" % template_group)

            else:
                template = template.replace("{{%s}}" % template_group, tree_data[template_group])

        with open(tree_path, 'w') as f:
            f.write(template)

    return tree_path


#
# Ensure that a given path is absolute.
#
def make_path_absolute(dir, base_dir):
    if base_dir not in dir:
        dir = os.path.join(base_dir, dir)

    return dir


#
# Initialise a database in the current folder.
#
def initialise_db(db_name):
    from expam.run import JSONConfig

    path = os.path.join(CWD, db_name)
    conf_path = os.path.join(path, "expam.conf")

    database_path = os.path.join(path, "database")
    phylogeny_path = os.path.join(path, "phylogeny")
    logs_path = os.path.join(path, "logs")
    results_path = os.path.join(path, "results")

    # Create the required folders.
    try:
        os.makedirs(phylogeny_path)

    except FileExistsError:
        # Check if the main folder already exists.
        if os.path.exists(conf_path):
            try:
                _conf = JSONConfig(url=conf_path)
                die("Database already exists!")

            except json.decoder.JSONDecodeError:
                die("Invalid name: non-database folder exists!")

        else:
            # Database has been initialised, but missing configuration file.
            # Get user confirmation this is indeed a database directory.
            is_db = get_user_confirmation("Database folder appears damaged. Restart it? (y/n)")
            if not is_db:
                die("Aborting...")

    except OSError:
        raise Exception("Failed to make directories %s!" % phylogeny_path)

    else:
        print("Successfully created directories %s!" % phylogeny_path)

    # Initialise the configuration file.
    conf = JSONConfig(name=db_name)
    conf.save(conf_path)
    print("Fresh database configuration generated!")

    # Create folder for logs.
    while True:
        try:
            os.mkdir(logs_path)
            print("Logs path %s created!" % logs_path)

            break

        except FileExistsError:
            print("Logs path %s already exists. Cleaning..." % logs_path)

            shutil.rmtree(logs_path)
            print("Previous logs cleared!")

    # Create database and results path.
    for db_path in (database_path, results_path):
        try:
            os.mkdir(db_path)

            print("Made path %s." % db_path)

        except FileExistsError:
            continue


#
# Modify sequences in database configuration file.
#
def _modify(directory, conf_path, first_n=None, group=None, add=True):
    from expam.run import JSONConfig
    from expam.sequences import ls

    # Ensure absolute path.
    base_url = make_path_absolute(
        dir=directory,
        base_dir=CWD
    )

    # Get list of sequences in base_url.
    seq_names = ls(base_url)

    if (first_n is None) or (first_n > len(seq_names)):
        first_n = len(seq_names)

    seq_dirs = [
        os.path.join(base_url, seq_names[i])
        if seq_names[i] != base_url else seq_names[i]
        for i in range(first_n)
    ]

    # Append directory to database configuration.
    conf = JSONConfig(url=conf_path)

    if add:
        conf.add_sequence(seq_dirs, group=group)
        print("Added %d files from %s." % (len(seq_dirs), base_url))
    else:
        conf.remove_sequence(seq_dirs, group=group)
        print("Removed %d files from %s." % (len(seq_dirs), base_url))

    conf.save()


def add_sequences(directory, conf_path, first_n=None, group=None):
    _modify(directory, conf_path, first_n, group, add=True)


def remove_sequences(directory, conf_path, first_n=None, group=None):
    _modify(directory, conf_path, first_n, group, add=False)


#
# Start building database with current configuration.
#
def build_database(db_path, plot):
    from expam.run import main as expam_main, JSONConfig, Timer

    conf_path = os.path.join(db_path, "expam.conf")
    conf = JSONConfig(conf_path)

    # Check if a phylogeny has been provided.
    if not check_phylogeny(conf):
        do_mashtree(conf, phy_dir=os.path.join(db_path, 'phylogeny'), group=None)

    # Read configuration file.
    k, n, phylogeny, genome_paths, pile = load_configuration_file(db_path)

    # Ensure correct phylogeny path.
    phylogeny = make_path_absolute(
        dir=phylogeny,
        base_dir=db_path
    )

    clear_logs(os.path.join(db_path, "logs"))

    # Set building database.
    try:
        with Timer() as t:
            expam_main(
                out_dir=db_path,
                genome_paths=genome_paths,
                phylogeny_path=phylogeny,
                k=k,
                n=n - 1,  # Account for main process.
                pile_size=pile,
                plot=plot
            )

        print("expam: " + str(t))

    finally:
        print(f'\nPID - {os.getpid()} dying...')


#
# Delete old log files.
#
def clear_logs(log_path):
    print("Clearing old log files...")

    try:
        shutil.rmtree(log_path)

    except FileNotFoundError:
        pass

    try:
        os.mkdir(log_path)

    except OSError:
        die("Can't make log path %s!" % log_path)


#
# Importable main.
#
def main():
    # Set process creation method.
    if platform.system() != "Windows":
        mp.set_start_method("fork")

    # Set command line arguments.
    parser = ArgumentParser(description="  expam CLI\n--------------\n", formatter_class=RawTextHelpFormatter)
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
    parser.add_argument("-n", "--n_processes", dest="n",
                        help="Number of CPUs to use for processing.",
                        metavar="[n (int)]")
    parser.add_argument("-s", "--sketch", dest="sketch",
                        help="Sketch size for mash.",
                        metavar="[sketch size (int)]")
    parser.add_argument("-p", "--phylogeny", dest="phylogeny",
                        help="URL of Newick file containing phylogeny.",
                        metavar="[phylogeny URL]")
    parser.add_argument("-d", "--directory", dest="directory",
                        help="File URL, context depending on command supplied.",
                        metavar="[directory]")
    parser.add_argument("-l", "--length", dest="length",
                        help="Length of simulated reads.",
                        metavar="[read length]")
    parser.add_argument("-o", "--out", dest="out_url",
                        help="URL to save sequences.",
                        metavar="[out URL]")
    parser.add_argument("-y", "--pile", dest="pile",
                        help="Number of genomes to pile at a time (or inf).",
                        metavar="[pile size]")
    parser.add_argument("-e", "--error_rate", dest="error_rate",
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
    parser.add_argument("--keep_zeros", dest="keep_zeros", default=False, action="store_true",
                        help="Keep nodes of output where no reads have been assigned.")
    parser.add_argument("--ignore_names", dest="ignore_node_names", default=False, action="store_true")
    parser.add_argument("--group", dest="groups", action="append", nargs="+",
                        help="Space-separated list of sample files to be treated as a single group in phylotree.")
    parser.add_argument("--colour_list", dest="colour_list", nargs="+",
                        help="List of colours to use when plotting groups in phylotree.")
    parser.add_argument("--name", dest="results_name", default=None,
                        help="Name of results folder.")
    parser.add_argument("--circle_scale", dest="circle_scale", default=1.0,
                        help="Scale of circles that represent splits in phylotree.")
    parser.add_argument("--sourmash", dest="use_sourmash", default=False, action="store_true",
                        help="Use sourmash for distance estimation.")
    parser.add_argument("--rapidnj", dest="use_rapidnj", default=True, action="store_true",
                        help="Use RapidNJ for Neighbour-Joining algorithm.")
    parser.add_argument("--quicktree", dest="use_quicktree", default=False, action="store_true",
                        help="Use QuickTree for Neighbour-Joining algorithm.")
    parser.add_argument("--paired", dest="paired_end", default=False, action="store_true",
                        help="Treat reads as paired-end.")
    parser.add_argument("--alpha", dest="alpha", default=None,
                        help="Percentage requirement for classification subtrees (see Tutorials 1 & 2).")

    # Parse arguments.
    args = parser.parse_args()

    # Argument groups.
    runtime_args = args.command, args.db_name, args.k, args.n, args.phylogeny, args.alpha
    directory_args = args.directory, args.out_url, args.truth_dir
    param_args = args.length, args.pile, args.error_rate, args.first_n, args.sketch, args.paired_end
    summary_args = args.plot, args.cutoff, args.cpm, args.taxonomy, args.results_name
    plot_args = args.groups, args.phyla, args.keep_zeros, not args.ignore_node_names, args.colour_list, \
                args.circle_scale, args.rank
    tree_args = args.use_sourmash, args.use_rapidnj, args.use_quicktree

    command, db_name, k, n, phylogeny, alpha = runtime_args
    directory, out_url, truth_dir = directory_args
    length, pile_size, error_rate, first_n, sketch, paired_end = param_args
    plot, cutoff, cpm, taxonomy, results_name = summary_args
    groups, plot_phyla, keep_zeros, use_node_names, colour_list, circle_scale, at_rank = plot_args
    use_sourmash, use_rapidnj, use_quicktree = tree_args

    group = None if groups is None else groups[0][0]  # When referring to sequence groups, not plotting groups.

    # Check we have a database name.
    if db_name is None and command not in AUTO_COMMANDS:
        try:
            DB_DIR = os.environ["EXPAM_DB_DIR"]
            db_name = os.path.basename(DB_DIR)

        except KeyError:
            die("A database name is required!")
            return
    else:
        DB_DIR = make_path_absolute(str(db_name), CWD)

    # Check we have a command.
    if command is None:
        die("No command supplied!")

    # Check valid form of passed parameters.
    k, n, length, first_n, cutoff = check_int_params(k, n, length, first_n, cutoff)
    error_rate, cpm, circle_scale = check_float_params(error_rate, cpm, circle_scale)

    # Determine supposed configuration file directory.
    CONF_DIR = os.path.join(DB_DIR, "expam.conf")
    PHY_DIR = os.path.join(DB_DIR, "phylogeny")

    def is_hex(string):
        if string[0] != "#" or len(string) != 7:
            return False

        try:
            int(string[1:].upper(), 16)
        except TypeError:
            return False

        return True

    # Format groups.
    if groups is not None:
        groups = [v for v in groups if v]

        for i in range(len(groups)):
            if is_hex(groups[i][0]):
                groups[i] = (groups[i][0], tuple(groups[i][1:]))
            else:
                groups[i] = (None, tuple(groups[i]))

    # Check colour list.
    if colour_list is not None:
        colour_list = [hex for hex in colour_list if is_hex(hex)]

        if not colour_list:
            colour_list = None

    # If plotting phyla, check we have taxonomic information.
    if (not plot_phyla) and (at_rank is None):
        name_to_lineage = None

    else:
        try:
            from expam.classification import load_taxonomy_map

            name_to_lineage, name_to_rank = load_taxonomy_map(DB_DIR)
        except FileNotFoundError:
            die("Run command `download_taxonomy` first to collect taxa for your genomes!")

            return

    # Command loop.

    #
    # Create and set configuration, then build database all in one step.
    #   *** Assumes all reference sequences are located in the same folder.
    #
    if command == "quickrun":
        from expam.run import JSONConfig

        # Initialise database folders.
        initialise_db(db_name)

        # Set the requested database parameters.
        conf = JSONConfig(url=CONF_DIR)
        conf.set(k=k, n=n, phylogeny_path=phylogeny, pile=pile_size)
        conf.save()

        # Add all sequences in the provided folder.
        add_sequences(directory, CONF_DIR, first_n)

        # Start building database.
        build_database(DB_DIR, plot)

    #
    # Print setup for default environment variable.
    #
    elif command == "default_db":
        print("export EXPAM_DB_DIR=%s" % DB_DIR)


    #
    # Initialise a database.
    #
    elif command == "create":
        initialise_db(db_name)

    #
    # Build database with set parameters.
    #
    elif command == "build":
        build_database(DB_DIR, plot)

    #
    # Print current parameter values in configuration file.
    #
    elif command == "print":
        # Check database exists.
        if not os.path.exists(CONF_DIR):
            die("'%s' not a database!" % db_name)

        config = load_configuration_file(DB_DIR, return_config=True)
        print(config)

    #
    # Run database against some reads.
    #
    elif command == "run":
        from expam.classification import run_classifier

        logs_path = os.path.join(DB_DIR, 'logs')
        clear_logs(logs_path)

        # Read configuration file.
        config = load_configuration_file(DB_DIR, return_config=True)

        k, n, phylogeny_path = int(config["k"]), int(config["n"]), config["phylogeny_path"]
        keys_shape, values_shape = tuple(config["keys_shape"]), tuple(config["values_shape"])

        # Correct for trimmed keys_shape.
        if len(keys_shape) == 1:
            keys_shape = keys_shape + (1,)

        phylogeny_path = make_path_absolute(phylogeny_path, DB_DIR)

        # Check alpha parameter is defined.
        alpha = calculate_alpha(k) if alpha is None else float(alpha)

        # Run expam classification.
        run_classifier(
            reads_url=directory,
            out_dir=DB_DIR,
            k=k,
            n=n - 1,  # Account for main process.
            phylogeny_path=phylogeny_path,
            keys_shape=keys_shape,
            values_shape=values_shape,
            logging_dir=logs_path,
            taxonomy=taxonomy,
            cutoff=cutoff,
            groups=groups,
            cpm=cpm,
            use_node_names=use_node_names,
            phyla=plot_phyla,
            name_taxa=name_to_lineage,
            colour_list=colour_list,
            results_name=results_name,
            circle_scale=circle_scale,
            paired_end=paired_end,
            alpha=alpha
        )

    #
    # Download mapping from accession ids to taxonomy.
    #
    elif command == "download_taxonomy":
        from expam.classification import accession_to_taxonomy

        # Convert accession ids to taxonomic lineages.
        accession_to_taxonomy(DB_DIR)

    #
    # Convert phylogenetic results to taxonomic results.
    #
    elif command == "to_taxonomy":
        if results_name is None:
            die("Require results directory (--name)!")

        from expam.classification import TAXID_LINEAGE_MAP_NAME, PHY_RESULTS, TAX_RESULTS, \
            load_taxonomy_map, name_to_id, ClassificationResults

        all_results_dir = os.path.join(DB_DIR, "results")
        results_dir = os.path.join(all_results_dir, results_name) \
            if results_name[:len(all_results_dir)] != all_results_dir else results_name

        if not os.path.exists(results_dir):
            die("Could not find results %s!" % results_dir)

        map_url = os.path.join(DB_DIR, TAXID_LINEAGE_MAP_NAME)
        if not os.path.exists(map_url):
            die("Run command `download_taxonomy` first to collect taxa for your genomes!")

        config = load_configuration_file(DB_DIR, return_config=True)
        phylogeny_path = config["phylogeny_path"]
        phylogeny_path = make_path_absolute(phylogeny_path, DB_DIR)

        index, phylogenyIndex = name_to_id(phylogeny_path)
        name_to_lineage, taxon_to_rank = load_taxonomy_map(DB_DIR)

        # Load phylogenetic results.
        phy_results_url = os.path.join(results_dir, PHY_RESULTS)
        tax_results_path = os.path.join(results_dir, TAX_RESULTS)

        if not os.path.exists(tax_results_path):
            os.mkdir(tax_results_path)

        results = ClassificationResults(
            index=index,
            phylogeny_index=phylogenyIndex,
            in_dir=phy_results_url,
            out_dir=results_dir,
            groups=groups,
            keep_zeros=keep_zeros,
            cutoff=cutoff,
            cpm=cpm,
            use_node_names=use_node_names,
            phyla=plot_phyla,
            name_taxa=name_to_lineage,
            colour_list=colour_list,
            circle_scale=circle_scale
        )
        results.to_taxonomy(name_to_lineage, taxon_to_rank, tax_results_path)
        results.draw_results()

    #
    # Plot tree.
    #
    elif command == "draw_tree":
        from ete3 import AttrFace, faces, Tree, TreeStyle, NodeStyle
        from expam.classification import DISC_COLOURS, PHYLA_COLOURS

        config = load_configuration_file(DB_DIR, return_config=True)
        phylogeny_path = config["phylogeny_path"]

        with open(make_path_absolute(phylogeny_path, DB_DIR), "r") as f:
            newick_string = f.read().strip()

        tree = Tree(newick_string, format=1)
        rank_colours = {}

        def layout(node):
            nonlocal name_to_lineage, name_to_rank

            if node.is_leaf():
                if use_node_names:
                    faces.add_face_to_node(AttrFace("name", fsize=20), node, column=0, position="aligned")

                if at_rank is not None:
                    lineage = name_to_lineage[node.name]

                    for name in lineage:
                        rank = name_to_rank[name][1]

                        if rank == at_rank:
                            if name not in rank_colours:
                                try:
                                    rank_colours[name] = DISC_COLOURS.pop()
                                except IndexError:
                                    die("Too many members of that rank!")

                            node_style = NodeStyle()
                            node_style['bgcolor'] = rank_colours[name]
                            node.set_style(node_style)

                elif plot_phyla:
                    for phyla, phyla_colour in PHYLA_COLOURS:
                        lineage = name_to_lineage[node.name]

                        if phyla in lineage:
                            node_style = NodeStyle()
                            node_style['bgcolor'] = phyla_colour
                            node.set_style(node_style)

        # Print tree to file.
        ts = TreeStyle()
        ts.mode = "c"
        ts.show_leaf_name = False
        ts.layout_fn = layout
        ts.force_topology = True
        ts.allow_face_overlap = False
        ts.draw_guiding_lines = True
        ts.root_opening_factor = 1

        phy_path = "phylotree.pdf" if not out_url else os.path.join(out_url)
        tree.render(
            phy_path,
            tree_style=ts,
            w=4000,
            units="px"
        )
        print("Phylogenetic tree written to %s!" % phy_path)

    #
    # Calculate precision and recall.
    #
    elif command == "evaluate":
        if not os.path.exists(directory):
            die("Couldn't find results path %s!" % directory)

        if not os.path.exists(truth_dir):
            die("Couldn't find truth dataset %s!" % truth_dir)

        import numpy as np
        import pandas as pd

        from expam.classification import load_taxonomy_map, request_labels, RANK_ORDER

        # Performance statistics.
        x = ['phylum', 'class', 'order', 'family', 'genus', 'species']
        y_prec_lib = []
        y_prec_cons = []
        y_rec_lib = []
        y_rec_cons = []

        taxid_to_lineage, taxon_to_rank = load_taxonomy_map(DB_DIR, convert_to_name=False)

        # Classifier results.
        expam_df = pd.read_csv(directory, sep="\t", names=['taxid', 'c_perc', 'c_cumul', 'c_count',
                                                           's_perc', 's_cumul', 's_count', 'rank', 'lineage'])
        expam_df['t_cumul'] = expam_df['s_cumul'] + expam_df['c_cumul']

        # Ensure we have representation for each rank in x.
        lineages = list({
            lineage
            for taxid, lineage in taxid_to_lineage.items()
            if taxid in expam_df['taxid']
        })

        for lineage in lineages:
            current_ranks = [(i, RANK_ORDER.index(taxon_to_rank[child][1])) for i, child in enumerate(lineage)]
            current_rank_indices = [rank_index for ind, rank_index in current_ranks]
            missing_rank_indices = [RANK_ORDER.index(rank) for rank in x if rank not in current_rank_indices]

            for rank_index in missing_rank_indices:
                # Find the closest child.
                child_index = min([ind for ind, child_rank_index in current_ranks
                                   if child_rank_index > rank_index])

                # Add this closest child as the representative member for this missing rank.
                child_taxid = taxon_to_rank[lineage[child_index]][0]
                child_df_index = expam_df.index[expam_df['taxid'] == child_taxid]

                new_row = expam_df.loc[child_df_index].copy()
                new_row['rank'] = RANK_ORDER[rank_index]

                # Append row to dataframe.
                expam_df.append(new_row, ignore_index=True)

        # Truth set.
        truth_df = pd.read_csv(truth_dir, sep="\t", names=['taxid', 'counts', 'abundance', 'rank', 'scientific name'])
        truth_df['taxid'] = truth_df['taxid'].apply(str)

        # Collect lineages for each taxid.
        species_taxids = truth_df['taxid'].tolist()
        taxid_to_lineage, taxon_to_rank = request_labels("taxonomy", "xml", species_taxids)

        taxid_to_lineage = {
            taxid: lineage.split(",")
            for taxid, lineage in taxid_to_lineage
        }
        taxon_to_rank = {
            taxon: rank.split(",")
            for taxon, rank in taxon_to_rank.items()
        }

        def at_rank(lineage, rank):
            nonlocal taxon_to_rank, RANK_ORDER

            z = RANK_ORDER.index(rank)
            descendants = [
                taxon_to_rank[name]
                for name in lineage
                if RANK_ORDER.index(taxon_to_rank[name][1]) >= z
            ]

            return descendants[0][0]

        rank_taxids = {
            rank: list({at_rank(lineage, rank) for lineage in taxid_to_lineage.values()})
            for rank in x
        }

        summary = np.sum if not taxonomy else len

        for rank, taxa in rank_taxids.items():
            df_subset = expam_df[expam_df['taxid'].isin(taxa)]

            # Total number of reads classified at given rank.
            mask = (expam_df['rank'] == rank) & (expam_df['c_cumul'] > cutoff)
            mask_comb = (expam_df['rank'] == rank) & (expam_df['t_cumul'] > cutoff)

            n_classified = summary(expam_df[mask]['c_cumul'])
            n_classified_combined = summary(expam_df[mask_comb]['t_cumul'])

            n_total = n if not taxonomy else len(taxa)

            # Precision and recall for confidently classified.
            y_prec_cons.append(
                summary(
                    df_subset[df_subset['c_cumul'] > cutoff]['c_cumul']
                ) / n_classified)
            y_rec_cons.append(
                summary(
                    df_subset[df_subset['c_cumul'] > cutoff]['c_cumul']
                ) / n_total)

            # Precision and recall under liberal classification regime.
            y_prec_lib.append(
                summary(
                    df_subset[df_subset['t_cumul'] > cutoff]['t_cumul']
                ) / n_classified_combined)
            y_rec_lib.append(
                summary(
                    df_subset[df_subset['t_cumul'] > cutoff]['t_cumul']
                ) / n_total)

        fig, (prec_ax, rec_ax) = plt.subplots(1, 2, figsize=(20, 10))

        prec_ax.set_title("Precision")
        prec_ax.plot(x, y_prec_lib, label="liberal", color="green", marker='o')
        prec_ax.plot(x, y_prec_cons, label="conservative", color="blue", marker='o')
        prec_ax.set_ylim([-0.05, 1.05])
        prec_ax.grid()

        rec_ax.set_title("Recall")
        rec_ax.plot(x, y_rec_lib, label="liberal", color="green", marker='o')
        rec_ax.plot(x, y_rec_cons, label="conservative", color="blue", marker='o')
        rec_ax.set_ylim([-0.05, 1.05])
        rec_ax.grid()

        prec_ax.legend(loc="best")
        rec_ax.legend(loc="best")

        # Save figure to disk.
        fig.savefig(os.path.join(out_url, "prec_recall.png"))

        # Save raw counts to disk.
        df = pd.DataFrame(index=x)

        df['prec_lib'] = y_prec_lib
        df['prec_cons'] = y_prec_cons
        df['rec_lib'] = y_rec_lib
        df['rec_cons'] = y_rec_cons
        df.to_csv(os.path.join(out_url, "prec_recall.csv"), sep="\t")

    #
    # Add sequence to configuration file.
    #
    elif command == "add":
        add_sequences(directory, CONF_DIR, group=group, first_n=first_n)

    #
    # Remove sequence from configuration file.
    #
    elif command == "remove":
        remove_sequences(directory, CONF_DIR, group=group, first_n=first_n)

    #
    # Set parameters in configuration file.
    #
    elif command == "set":
        from expam.run import JSONConfig

        conf = JSONConfig(url=CONF_DIR)

        if group is None:
            conf.set(k=k, n=n, phylogeny_path=phylogeny, pile=pile_size, s=sketch)
        else:
            conf.group_set(group_name=group, k=k, s=sketch)

        conf.save()

    #
    # Sketch sequences.
    #
    elif command == "sketch":
        from expam.run import JSONConfig

        conf = JSONConfig(url=CONF_DIR)
        do_sketches(conf, PHY_DIR, group, use_sourmash=use_sourmash)


    #
    # Compute distance tree.
    #
    elif command == "distance":
        from expam.run import JSONConfig

        conf = JSONConfig(url=CONF_DIR)
        do_distances(conf, PHY_DIR, group, use_sourmash=use_sourmash)

    #
    # Execute neighbour joining on tree.
    #
    elif command == "nj":
        from expam.run import JSONConfig

        conf = JSONConfig(url=CONF_DIR)
        do_trees(conf, PHY_DIR, group, use_quicktree=use_quicktree)

    #
    # Build phylogeny.
    #
    elif command == "tree":
        from expam.run import JSONConfig

        def argmax(*checks):
            if not any(checks):
                return 0

            return max(i for i, val in enumerate(checks + (True, )) if val is True)

        conf = JSONConfig(url=CONF_DIR)

        entry_points = (
            partial(do_sketches, conf=conf, phy_dir=PHY_DIR,
                    group=group, use_sourmash=use_sourmash),    # Sketch sequences.
            partial(do_distances, conf=conf, phy_dir=PHY_DIR,
                    group=group, use_sourmash=use_sourmash),    # Pairwise distances.
            partial(do_trees, conf=conf, phy_dir=PHY_DIR,
                    group=group, use_quicktree=use_quicktree)   # NJ trees.
        )
        entry_stage = argmax(
            check_sketches(conf, PHY_DIR, group, use_sourmash=use_sourmash),
            check_distances(conf, PHY_DIR, group),
            check_trees(conf, PHY_DIR, group)
        )

        for stage in range(entry_stage, 3):
            entry_points[stage]()

        tree_dir = finalise_tree(conf, PHY_DIR, db_name)
        conf.set(phylogeny_path=tree_dir)
        conf.save()

        print("Phylogeny set to %s." % tree_dir)

    #
    # Make mashtree for database.
    #
    elif command == "mashtree":
        from expam.run import JSONConfig

        conf = JSONConfig(url=CONF_DIR)
        do_mashtree(conf, PHY_DIR, group)

    #
    # Create reads from some reference sequence. Used for unit tests.
    #
    elif command == "make_reads":
        from expam.sequences import make_reads

        make_reads(
            in_url=directory,
            l=length,
            n=n,
            out_url=out_url,
            e=0.0 if error_rate is None else float(error_rate)
        )

    #
    # Posthumously employ cutoff.
    #
    elif command == "cutoff":
        from expam.sequences import ls

        if os.path.exists(results_name):
            results_dir = results_name
        else:
            results_dir = os.path.join(DB_DIR, "results", results_name)

        tax_results_dir = os.path.join(results_dir, "tax")
        raw_results_dir = os.path.join(tax_results_dir, "raw")

        if not os.path.exists(results_dir):
            die("Could not find results %s!" % results_dir)

        file_names = [
            file_dir
            for file_dir in ls(tax_results_dir)
        ]

        for summary_file in file_names:
            print("Employing cutoff on sample %s..." % summary_file)

            with open(os.path.join(tax_results_dir, summary_file), "r") as f:
                data = [line.strip().split('\t') for line in f]

            total_classifications = sum(int(line[3]) for line in data[1:])
            total_splits = sum(int(line[6]) for line in data[1:])
            total_reads = total_classifications + total_splits

            if cpm is not None:
                cutoff = cpm * (total_reads / 1e6)
            elif cutoff is None:
                raise Exception("Require a cutoff value!")

            def valid_line(line):
                if int(line[2]) >= cutoff or int(line[5]) >= cutoff:
                    return True

            valid_taxids = [line[0] for line in data[1:] if valid_line(line)]
            valid_data = [line for i, line in enumerate(data) if i == 0 or valid_line(line)]

            # Overwrite summary.
            with open(os.path.join(tax_results_dir, summary_file), "w") as f:
                f.write("\n".join("\t".join(line) for line in valid_data))

            # Overwrite raw output.
            with open(os.path.join(raw_results_dir, summary_file), "r") as f:
                valid_data = []

                for line in f:
                    current_data = line.strip().split('\t')
                    if current_data[2] in valid_taxids:
                        valid_data.append('\t'.join(current_data))

            with open(os.path.join(raw_results_dir, summary_file), "w") as f:
                f.write('\n'.join(valid_data))

    #
    # Generate evenly-balanced phylogeny for sequences contained
    #   in directory.
    #
    elif command == "fake_phylogeny":
        from expam.run import JSONConfig
        from expam.sequences import format_name
        from expam.tree import simulate_balanced_phylogeny

        PHY_NAME = "expam_outtree.nwk"
        PHY_DIR = os.path.join("phylogeny", PHY_NAME)

        # Load configuration file.
        config = load_configuration_file(DB_DIR, return_config=True)

        # Get list of sequences.
        seq_names = [
            format_name(os.path.basename(url))
            for url in config.get_paths()
        ]

        if not seq_names:
            die("No sequences added to the database!")

        # Get valid save_directory.
        save_dir = os.path.join(DB_DIR, PHY_DIR)

        # Compute balanced phylogeny.
        newick_str = simulate_balanced_phylogeny(seq_names)

        with open(save_dir, "w") as f:
            f.write(newick_str)

        print("Newick file generated at %s." % save_dir)

        # Add phylogeny to expam configuration file.
        conf = JSONConfig(url=CONF_DIR)
        conf.set(phylogeny_path=PHY_DIR)
        conf.save()

        print("Database phylogeny set to %s." % save_dir)

    #
    # Count k-mers in sequence file.
    #
    elif command == "count":
        import numpy as np
        import pandas as pd

        from expam.c.extract import extract_kmers_from_read as get_kmers
        from expam.sequences import format_name
        from expam.stores import SequenceStore

        store = SequenceStore(k=k, out_dir="")

        seq_name = format_name(directory)
        store.add_sequence(directory)

        kmers = get_kmers(store[seq_name], k)
        values, counts = np.unique(kmers, return_counts=True)

        df = pd.DataFrame({"kmers": values, "counts": counts})
        df.to_csv(out_url, sep="\t", index=False, header=False)

    #
    # Map Jellyfish counts to expam.
    #
    elif command == "from_jellyfish":
        import numpy as np
        import pandas as pd

        from expam.c.extract import map_kmers
        from expam.run import calculate_n_chunks, KEYS_DATA_TYPE

        df = pd.read_csv(directory, sep=" ", names=["kmers", "counts"])
        seq = ("N".join(df['kmers'])).encode("utf8")

        kmers = np.zeros((len(df['kmers']), calculate_n_chunks(k) + 1), dtype=KEYS_DATA_TYPE)
        kmers[:, 0] = np.arange(len(df['kmers']), dtype=KEYS_DATA_TYPE) * (k + 1)
        map_kmers(seq, k, kmers)

        for i in range(calculate_n_chunks(k)):
            df['kmers_%d' % i] = kmers[:, i + 1]

        df.drop(columns=['kmers'], inplace=True)
        df.to_csv(out_url, sep="\t", index=False, header=False,
                  columns=['kmers_%d' % i for i in range(calculate_n_chunks(k))] + ['counts'])

    #
    # Plot time taken for tasks.
    #
    elif command == "plot":
        from expam.run import plot_data

        logs_dir = os.path.join(DB_DIR, "logs")
        _, t0 = get_t0(logs_dir)

        plot_data(logs_dir, t0)

    #
    #   Plot memory used by each process over time.
    #
    elif command == "plot_memory":
        from expam.sequences import ls

        def _get_mem(data):
            start, end = data.find("(") + 1, data.find(")")

            if start == -1 or end == -1:
                raise ValueError("Invalid data line %s!" % data)

            rss_string, vsz_string = data[start:end].split(", ")

            for point in (rss_string, vsz_string):
                mem_type, mem_used = point.split(" ")

                yield float(mem_used.replace("Gb", ""))

        def _r_yield_dict(obj, current_key=""):
            if isinstance(obj, dict):
                for key, value in obj.items():
                    yield from _r_yield_dict(value, "%s_%s" % (current_key, str(key)))

            else:
                yield current_key, obj

        job_types = {
            "_extract.log": "extract",
            "_union.log": "union",
        }
        plot_colours = {
            "extract": {
                "rss": "blue",
                "vsz": "black",
            },
            "union": {
                "rss": "green",
                "vsz": "purple"
            }
        }
        job_labels = {
            "extract": "extract",
            "union": "union",
        }

        logs_dir = os.path.join(DB_DIR, "logs")

        log_urls = [
            os.path.join(logs_dir, log_name)
            for log_name in ls(logs_dir, ext=".log")
            if "_main" not in log_name
        ]

        _, t0 = get_t0(logs_dir)

        # Collect total memory usage as well.
        cumulative = {
            "union": {
                "rss": 0,
                "vsz": 0,
            },
            "extract": {
                "rss": 0,
                "vsz": 0,
            }
        }

        fig, (memory_time_ax, cumulative_memory_ax) = plt.subplots(1, 2, figsize=(20, 10))

        for log_url in log_urls:

            with open(log_url, "r") as f:
                log_data = f.readlines()

            t, rss, vsz = [], [], []

            for line in log_data:
                if "of RAM" in line:
                    rss_mem, vsz_mem = _get_mem(line)
                    time_point = (get_time(line) - t0).total_seconds()

                    rss.append(rss_mem)
                    vsz.append(vsz_mem)
                    t.append(time_point / 3600)  # Time in hours.

            for mem_type, mem_data in (("rss", rss), ("vsz", vsz)):
                for k, v in job_types.items():
                    if k in log_url:
                        job_type = v

                        break
                else:
                    raise TypeError("Unknown process log type %s!" % os.path.basename(log_url))

                plot_label = job_labels[job_type]
                # Don't allow duplicate labels.
                job_labels[job_type] = ""

                memory_time_ax.plot(
                    t,
                    mem_data,
                    label=plot_label,
                    marker=".",
                    color=plot_colours[job_type][mem_type]
                )

                if mem_data:
                    # Add to cumulative memory usage.
                    cumulative[job_type][mem_type] += mem_data[-1]

        # Memory-over-time axis labels.
        memory_time_ax.set(
            xlabel="Time since start [hrs]",
            ylabel="Memory used [Gb]",
            title="Memory used by child workers"
        )
        memory_time_ax.legend(loc="best")

        # Bar chart of cumulative memory usage.
        labels = ["union", "extract"]
        rss = [cumulative["union"]["rss"], cumulative["extract"]["rss"]]
        vsz = [cumulative["union"]["vsz"], cumulative["extract"]["vsz"]]

        # Set axis positions.
        x = list(range(len(labels)))
        width = 0.35

        rss_x = list(v - width / 2 for v in x)
        vsz_x = list(v + width / 2 for v in x)

        # Plot bar charts.
        rss_bars = cumulative_memory_ax.bar(rss_x, rss, width, label="RSS")
        vsz_bars = cumulative_memory_ax.bar(vsz_x, vsz, width, label="VSZ")

        # Cumulative memory usage axis labels.
        cumulative_memory_ax.set(
            ylabel="Total memory used [Gb]",
            xticks=x,
            xticklabels=labels,
            title="Total memory used by workers",
        )
        cumulative_memory_ax.legend(loc="best")

        # Save figure to disk.
        fig.savefig(os.path.join(logs_dir, "memory_usage.png"))

    #
    #   Count how many genomes were processed in a given run.
    #
    elif command == "how_many":
        from expam.sequences import ls

        def _get_log_lines(log_dir):
            with open(log_dir, "r") as f:
                lines = f.readlines()

            return lines

        logs_dir = os.path.join(DB_DIR, "logs")
        log_files = [
            os.path.join(logs_dir, log_name)
            for log_name in ls(logs_dir, ext=".log")
        ]

        # First get extractions done.
        genomes_extracted = 0

        for log_file in log_files:
            if "_extract.log" in log_file:
                log_lines = _get_log_lines(log_file)

                for line in log_lines:
                    if "kmers sent!" in line:
                        genomes_extracted += 1

                    # Only count those during extraction phase.
                    elif "now at phase map!" in line:
                        break

        print("Genome k-mers extracted:\t%d." % genomes_extracted)

        # Now count how many unions done.
        genomes_unioned = []

        for log_file in log_files:
            if "_union.log" in log_file:
                log_lines = _get_log_lines(log_file)

                log_unions = 0
                for line in log_lines:
                    if "Importing" in line:
                        log_unions += 1

                genomes_unioned.append(log_unions)

        print("Genomes undergone union:\t%d." % min(genomes_unioned))
        print("\t Individual unions ~ {%s}" % " ".join(str(union) for union in genomes_unioned))

    else:
        die("Invalid command: " + str(command) + ".")


if __name__ == "__main__":
    main()
