import os
from typing import Set
import matplotlib.pyplot as plt
from expam.classify import ResultsPathConfig
from expam.classify.config import make_results_config, validate_classification_results, validate_results_configuration
from expam.classify.taxonomy import TaxonomyNCBI
from expam.cli.main import CommandGroup, ExpamOptions
from expam.database import FileLocationConfig
from expam.database.config import JSONConfig, make_database_config, validate_database_file_configuration
from expam.sequences import format_name
from expam.tree.simulate import simulate_balanced_phylogeny
from expam.utils import die, ls


class UtilsCommand(CommandGroup):
    commands: Set[str] = {
        'download_taxonomy', 'cutoff', 'fake_phylogeny', 'plot_memory'
    }

    def __init__(
        self, config: FileLocationConfig, out_dir: str, cutoff: int, cpm: float
    ) -> None:
        super().__init__()

        self.config: FileLocationConfig = config
        self.out_dir: str = out_dir
        self.results_config: ResultsPathConfig = None if out_dir is None else make_results_config(out_dir)

        self.json_conf = None

        self.count_cutoff = cutoff
        self.cpm = cpm

    @classmethod
    def take_args(cls, args: ExpamOptions) -> dict:
        cutoff = cls.parse_ints(args.cutoff)
        cpm = cls.parse_floats(args.cpm)

        return {
            'config': make_database_config(args.db_name),
            'out_dir': args.out_url,
            'cutoff': cutoff,
            'cpm': cpm
        }

    def get_conf(self):
        if self.json_conf is None:
            self.check_database_exists()
            self.json_conf = JSONConfig(self.config.conf)
        
        return self.json_conf

    def check_database_exists(self):
        validate_database_file_configuration(self.config)

    """
    Download taxonomy command
    =========================
    
    """
    def download_taxonomy(self):
        self.check_database_exists()

        tax_obj: TaxonomyNCBI = TaxonomyNCBI(self.config)
        tax_obj.accession_to_taxonomy()

    """
    Employ cutoff on taxonomic classification output
    ================================================
    
    """
    def cutoff(self):
        if self.out_dir is None:
            die("Must supply -o/--out!")

        validate_results_configuration(self.results_config, check_taxonomy=True)

        if os.path.exists(self.out_dir):
            validate_classification_results(self.out_dir)
        else:
            raise ValueError("Results path %s doesn't exist!")

        file_names = [
            file_dir
            for file_dir in ls(self.results_config.tax)
            if file_dir != self.results_config.tax_classified and file_dir != self.results_config.tax_split
        ]

        for summary_file in file_names:
            print("Employing cutoff on sample %s..." % summary_file)

            with open(summary_file, "r") as f:
                data = [line.strip().split('\t') for line in f]

            total_classifications = sum(int(line[3]) for line in data[1:])
            total_splits = sum(int(line[6]) for line in data[1:])
            total_reads = total_classifications + total_splits

            if self.cpm is not None:
                cutoff = self.cpm * (total_reads / 1e6)
            elif cutoff is None:
                raise Exception("Require a cutoff value!")

            def valid_line(line):
                if int(line[2]) >= cutoff or int(line[5]) >= cutoff:
                    return True

            valid_taxids = [line[0] for line in data[1:] if valid_line(line)]
            valid_data = [line for i, line in enumerate(data) if i == 0 or valid_line(line)]

            # Overwrite summary.
            with open(summary_file, "w") as f:
                f.write("\n".join("\t".join(line) for line in valid_data))

            # Overwrite raw output.
            raw_sample_output: str = os.path.join(self.results_config.tax_raw, os.path.basename(summary_file))
            with open(raw_sample_output, "r") as f:
                valid_data = []

                for line in f:
                    current_data = line.strip().split('\t')
                    if current_data[2] in valid_taxids:
                        valid_data.append('\t'.join(current_data))

            with open(raw_sample_output, "w") as f:
                f.write('\n'.join(valid_data))

    """
    Generate a fake phylogeny
    =========================
    
    """
    def fake_phylogeny(self):
        conf: JSONConfig = self.get_conf()

        # Get list of sequences.
        seq_names = [
            format_name(os.path.basename(url))
            for url in conf.get_paths()
        ]

        if not seq_names:
            die("No sequences added to the database!")

        # Compute balanced phylogeny.
        fake_phy_path = os.path.join(self.config.phylogeny, 'expam_outtree.nwk')
        newick_str = simulate_balanced_phylogeny(seq_names)

        with open(fake_phy_path, "w") as f:
            f.write(newick_str)

        print("Newick file generated at %s." % fake_phy_path)

        # Add phylogeny to expam configuration file.
        conf.set(phylogeny_path=fake_phy_path)
        conf.save()

        print("Database phylogeny set to %s." % fake_phy_path)

    """
    Plot expam database build memory usage
    ======================================
    
    """
    def plot_memory(self):
        conf: JSONConfig = self.get_conf()

        def _get_mem(data: str):
            start, end = data.find("(") + 1, data.find(")")

            if start == -1 or end == -1:
                raise ValueError("Invalid data line %s!" % data)

            rss_string, vsz_string = data[start:end].split(", ")

            for point in (rss_string, vsz_string):
                _, mem_used = point.split(" ")

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

        log_urls = [
            os.path.join(self.config.logs, log_name)
            for log_name in ls(self.config.logs, ext=".log")
            if "_main" not in log_name
        ]

        _, t0 = self.get_t0(self.config.logs)

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
                    time_point = (self.get_time(line) - t0).total_seconds()

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
        fig.savefig(os.path.join(self.config.logs, "memory_usage.png"))

    