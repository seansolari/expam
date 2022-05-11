import os
from typing import List, Set, Tuple
from expam.classify import ResultsPathConfig
from expam.classify.classify import ClassificationResults, name_to_id, run_classifier
from expam.classify.config import make_results_config, validate_results_configuration
from expam.classify.taxonomy import TaxonomyNCBI
from expam.cli.main import CommandGroup, ExpamOptions, clear_logs
from expam.database import FileLocationConfig
from expam.database.config import JSONConfig, make_database_config, validate_database_file_configuration
from expam.utils import die, is_hex, make_path_absolute


class ClassifyCommand(CommandGroup):
    commands: Set[str] = {
        'classify', 'to_taxonomy'
    }

    def __init__(
        self, config: FileLocationConfig,
        files: List[str], out_dir: str,
        convert_to_taxonomy: bool, cutoff: int, cpm: float, groups: List[Tuple[str]],
        use_node_names: bool, keep_zeros: bool, plot_phyla: bool,
        colour_list: List[str], paired_end: bool, alpha: float,
        log_scores: bool, itol_mode: bool
    ) -> None:
        super().__init__()

        self.config: FileLocationConfig = config
        self.results_config: ResultsPathConfig = make_results_config(out_dir)

        self.files = files
        self.out_dir = out_dir

        self.convert_to_taxonomy = convert_to_taxonomy
        self.cutoff = cutoff
        self.cpm = cpm
        self.groups = groups

        self.use_node_names = use_node_names
        self.keep_zeros = keep_zeros
        self.plot_phyla = plot_phyla
        self.colour_list = colour_list

        self.paired_end = paired_end
        self.alpha = alpha
        self.log_scores = log_scores
        self.itol_mode = itol_mode

    @classmethod
    def take_args(cls: CommandGroup, args: ExpamOptions) -> dict:
        cutoff = cls.parse_ints(args.cutoff)
        cpm, alpha = cls.parse_floats(args.cpm, args.alpha)

        if args.out_url is None:
            die("Must supply -o/--out.")

        # Format groups.
        if args.groups is not None:
            groups = [v for v in args.groups if v]

            for i in range(len(groups)):
                if is_hex(groups[i][0]):
                    groups[i] = (groups[i][0], tuple(groups[i][1:]))
                else:
                    groups[i] = (None, tuple(groups[i]))
        else:
            groups = None

        # Check colour list.
        if args.colour_list is not None:
            colour_list = [hex for hex in args.colour_list if is_hex(hex)]

            if not colour_list:
                colour_list = None
        else:
            colour_list = None

        return {
            'config': make_database_config(args.db_name),
            'files': args.directory,
            'out_dir': args.out_url,
            'convert_to_taxonomy': args.taxonomy,
            'cutoff': cutoff,
            'cpm': cpm,
            'groups': groups,
            'use_node_names': not args.ignore_names,
            'keep_zeros': args.keep_zeros,
            'plot_phyla': args.phyla,
            'colour_list': colour_list,
            'paired_end': args.paired_end,
            'alpha': alpha,
            'log_scores': args.log_scores,
            'itol_mode': args.itol_mode
        }

    def check_database_exists(self):
        validate_database_file_configuration(self.config)

        if not os.path.exists(self.config.database_file):
            die("Database has not been built! Not found at %s." % self.config.database_file)

    """
    Run command
    ===========
    
    """
    def classify(self):
        self.check_database_exists()
        clear_logs(self.config.logs)

        # Read configuration file.
        config: JSONConfig = JSONConfig(self.config.conf)
        k, n, phylogeny_path, _, _ = config.get_build_params()
        keys_shape, values_shape = tuple(config["keys_shape"]), tuple(config["values_shape"])

        # Correct for trimmed keys_shape.
        if len(keys_shape) == 1:
            keys_shape = keys_shape + (1,)

        try:
            tax_obj = TaxonomyNCBI(self.config)
            name_to_lineage, _ = tax_obj.load_taxonomy_map(self.config)
        except OSError:
            if self.convert_to_taxonomy:
                die("First run `download_taxonomy` to collect associated taxonomy data.")
            else:
                name_to_lineage = None

        # Run expam classification.
        run_classifier(
            read_paths=self.files,
            out_dir=self.out_dir,
            db_dir=self.config.base,
            k=k,
            n=n - 1,  # Account for main process.
            phylogeny_path=phylogeny_path,
            keys_shape=keys_shape,
            values_shape=values_shape,
            logging_dir=self.config.logs,
            taxonomy=self.convert_to_taxonomy,
            cutoff=self.cutoff,
            groups=self.groups,
            cpm=self.cpm,
            use_node_names=self.use_node_names,
            phyla=self.plot_phyla,
            name_taxa=name_to_lineage,
            colour_list=self.colour_list,
            paired_end=self.paired_end,
            alpha=self.alpha,
            log_scores=self.log_scores,
            itol_mode=self.itol_mode
        )
    
    """
    Convert phylogenetic results to taxonomy
    ========================================
    
    """
    def to_taxonomy(self):
        self.check_database_exists()

        if self.out_dir is None:
            die("Require output directory (-o, --out_dir)!")
        else:
            validate_results_configuration(self.results_config)

        if not os.path.exists(self.config.taxid_lineage):
            die("Run command `download_taxonomy` first to collect taxa for your genomes!")

        config = JSONConfig(self.config.conf)
        phylogeny_path = make_path_absolute(config["phylogeny_path"], self.config.conf)

        index, phylogenyIndex = name_to_id(phylogeny_path)

        tax_obj = TaxonomyNCBI(self.config)
        name_to_lineage, taxon_to_rank = tax_obj.load_taxonomy_map(self.config)

        if not os.path.exists(self.results_config.tax):
            os.mkdir(self.results_config.tax)

        results = ClassificationResults(
            index=index,
            phylogeny_index=phylogenyIndex,
            results_config=self.results_config,
            groups=self.groups,
            keep_zeros=self.keep_zeros,
            cutoff=self.cutoff,
            cpm=self.cpm,
            use_node_names=self.use_node_names,
            phyla=self.plot_phyla,
            name_taxa=name_to_lineage,
            colour_list=self.colour_list,
            log_scores=self.log_scores
        )
        results.to_taxonomy(name_to_lineage, taxon_to_rank)
    