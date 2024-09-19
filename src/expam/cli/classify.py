# version: 2
import logging
import os
from typing import List, Set, Tuple
from expam.classify.classify import TaxonomicResults, run_classifier
from expam.classify.config import create_tax_results, make_results_config, validate_results_configuration
from expam.classify.taxonomy import TaxonomyInterface
from expam.cli.main import CommandGroup, ExpamOptions, clear_logs
from expam.database import FileLocationConfig
from expam.database.config import JSONConfig, make_database_config, validate_database_file_configuration
from expam.tree.tree import Index
from expam.utils import die, is_hex
import expam.logger

class ClassifyCommand(CommandGroup):
    commands: Set[str] = {
        'classify', 'to_taxonomy'
    }

    def __init__(
        self, config: FileLocationConfig,
        files: List[str], out_dir: str,
        convert_to_taxonomy: bool, cpm: float, groups: List[Tuple[str]],
        use_node_names: bool, keep_zeros: bool,
        colour_list: List[str], paired_end: bool, alpha: float,
        log_scores: bool, itol_mode: bool, flat_colour: bool, debug: bool
    ) -> None:
        super().__init__()

        self.config = config
        self.results_config = make_results_config(out_dir)

        self.files = files
        self.out_dir = out_dir

        self.convert_to_taxonomy = convert_to_taxonomy
        self.cpm = cpm
        self.groups = groups

        self.use_node_names = use_node_names
        self.keep_zeros = keep_zeros
        self.colour_list = colour_list

        self.paired_end = paired_end
        self.alpha = alpha
        self.log_scores = log_scores
        self.itol_mode = itol_mode
        self.flat_colour = flat_colour

        if debug:
            expam.logger.current_logging_level = logging.DEBUG

    @classmethod
    def take_args(cls: CommandGroup, args: ExpamOptions) -> dict:
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
            'cpm': cpm,
            'groups': groups,
            'use_node_names': not args.ignore_names,
            'keep_zeros': args.keep_zeros,
            'colour_list': colour_list,
            'paired_end': args.paired_end,
            'alpha': alpha,
            'log_scores': args.log_scores,
            'itol_mode': args.itol_mode,
            'flat_colour': args.flat_colour,
            'debug': args.debug
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
        config = JSONConfig(self.config.conf)
        k, n, phylogeny_path, _, _ = config.get_build_params()
        keys_shape, values_shape = tuple(config["keys_shape"]), tuple(config["values_shape"])

        # Correct for trimmed keys_shape.
        if len(keys_shape) == 1:
            keys_shape = keys_shape + (1,)

        # Run expam classification.
        run_classifier(self.files, self.out_dir, self.config.base, k, n-1, phylogeny_path,
                       keys_shape, values_shape, self.config.logs,
                       taxonomy=self.convert_to_taxonomy, groups=self.groups,
                       keep_zeros=self.keep_zeros, cpm=self.cpm, use_node_names=self.use_node_names,
                       colour_list=self.colour_list, paired_end=self.paired_end, alpha=self.alpha,
                       log_scores=self.log_scores, itol_mode=self.itol_mode, flat_colour=self.flat_colour)
    
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
        # load phylogenetic tree
        config = JSONConfig(self.config.conf)
        phylogeny_path = config.get_phylogeny_path()
        _, tree = Index.load_newick(phylogeny_path)
        # create phylogenetic results
        create_tax_results(self.results_config)
        my_ncbi = TaxonomyInterface(self.config, tree=tree)
        tax = TaxonomicResults.from_tables(self.results_config.phy_classified,
                                           self.results_config.phy_split,
                                           tree,
                                           my_ncbi)
        tax.summarise(self.results_config.tax,
                      cpm=self.cpm)
        tax.convert_raw(self.results_config.phy_raw, self.results_config.tax_raw, tree)
