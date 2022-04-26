import os
import re
import shutil
import subprocess
from typing import List, Set, Tuple, Union
from expam.classify import ResultsPathConfig
from expam.classify.classify import ClassificationResults, name_to_id
from expam.classify.config import make_results_config, validate_results_configuration
from expam.classify.taxonomy import TaxonomyNCBI
from expam.cli.main import CommandGroup, ExpamOptions
from expam.database import FileLocationConfig
from expam.database.config import JSONConfig, make_database_config, validate_database_file_configuration
from expam.tree import PHYLA_COLOURS
from expam.tree.tree import RandomColour
from expam.utils import die, is_hex, ls, make_path_absolute


class TreeCommand(CommandGroup):
    commands: Set[str] = {
        'phylotree', 'draw_tree', 'tree', 
        'mashtree', 'sketch', 'distance', 'nj'
    }

    def __init__(
        self, config: FileLocationConfig,
        out_dir: str, cutoff: int, cpm: float, groups: List[Tuple[str]],
        use_node_names: bool, keep_zeros: bool, plot_phyla: bool,
        colour_list: List[str], log_scores: bool, itol_mode: bool,
        at_rank: str, use_sourmash: bool, use_quicktree: bool
    ) -> None:
        super().__init__()

        self.config: FileLocationConfig = config
        self.results_config: ResultsPathConfig = None if out_dir is None else make_results_config(out_dir)
        self.json_conf = None
        self.out_dir = out_dir

        self.cutoff = cutoff
        self.cpm = cpm
        self.group = None if groups is None else groups[0][0]
        self.groups = groups

        self.use_node_names = use_node_names
        self.keep_zeros = keep_zeros
        self.plot_phyla = plot_phyla
        self.colour_list = colour_list

        self.log_scores = log_scores
        self.itol_mode = itol_mode

        self.at_rank = at_rank

        self.use_sourmash = use_sourmash
        self.use_quicktree = use_quicktree


    @classmethod
    def take_args(cls, args: ExpamOptions) -> dict:
        cutoff, n = cls.parse_ints(args.cutoff, args.n)
        cpm = cls.parse_floats(args.cpm)

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
            'out_dir': args.out_url,
            'cutoff': cutoff,
            'cpm': cpm,
            'groups': groups,
            'use_node_names': not args.ignore_names,
            'keep_zeros': args.keep_zeros,
            'plot_phyla': args.phyla,
            'colour_list': colour_list,
            'log_scores': args.log_scores,
            'itol_mode': args.itol_mode,
            'at_rank': args.rank,
            'use_sourmash': args.use_sourmash,
            'use_quicktree': args.use_quicktree
        }

    def get_conf(self):
        if self.json_conf is None:
            self.check_database_exists()
            self.json_conf = JSONConfig(self.config.conf)
        
        return self.json_conf

    def get_group_or_die(self, group_name):
        conf: JSONConfig = self.get_conf()
        k, s, sequences = conf.group_get(group_name)

        if k is None:
            die("Parameter `k` has not been set for group %s." % group_name)
        elif s is None:
            die("Parameter `s` has not been set for group %s." % group_name)

        return k, s, sequences

    def get_n_processes(self) -> int:
        conf: JSONConfig = self.get_conf()
        n_processes: Union[None, int] = conf.get_n_processes()

        if n_processes is None:
            return 1
        else:
            return int(n_processes)

    def check_database_exists(self):
        validate_database_file_configuration(self.config)

    def check_results_exist(self):
        if self.results_config is None:
            die("Results not specified!")

        if not validate_results_configuration(self.config):
            die("Invalid results configuration!")

    @staticmethod
    def _check_command(command):
        __has = shutil.which(command) is not None
        if not __has:
            die('Could not find local installation of %s!' % command)

    def check_mash(self):
        self._check_command('mash')

    def check_mashtree(self):
        self._check_command('mashtree')

    def check_rapidnj(self):
        self._check_command('rapidnj')

    def check_quicktree(self):
        self._check_command('quicktree')

    @staticmethod
    def shell(command, cwd=None):
        if cwd is None:
            cwd = os.getcwd()

        p = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            cwd=cwd
        )
        if p.returncode != 0:
            raise Exception("\n" + p.stderr.decode("utf-8"))

        return p.stdout.decode("utf-8")
    
    """
    Make phylotree of classification output
    =======================================
    
    """
    def phylotree(self):
        self.check_results_exist()
        config: JSONConfig = self.get_conf()

        phylogeny_path = make_path_absolute(config["phylogeny_path"], self.config.database)
        index, phylogenyIndex = name_to_id(phylogeny_path)

        try:
            tax_obj: TaxonomyNCBI = TaxonomyNCBI(self.config)
            name_to_lineage, _ = tax_obj.load_taxonomy_map()
        except FileNotFoundError:
            name_to_lineage = None

        # Load phylogenetic results.
        results = ClassificationResults(
            index=index,
            phylogeny_index=phylogenyIndex,
            in_dir=self.results_config.phy,
            out_dir=self.out_dir,
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
        results.draw_results(itol_mode=self.itol_mode)

    """
    Draw phylogenetic tree
    ======================
    
    """
    def draw_tree(self):
        try:
            from ete3 import AttrFace, faces, Tree, TreeStyle, NodeStyle
        except ImportError:
            die("First install ete3 tools...`pip install ete3`.")
        
        config: JSONConfig = self.get_conf()
        try:
            phylogeny_path = make_path_absolute(config["phylogeny_path"], self.config.database)
        except TypeError:
            die("Phylogeny not set!")

        try:
            with open(phylogeny_path, "r") as f:
                newick_string = f.read().strip()
        except OSError:
            die("Can't find Newick file %s." % phylogeny_path)

        tree = Tree(newick_string, format=1)
        rank_colours = {}

        tax_obj: TaxonomyNCBI = TaxonomyNCBI(self.config)
        name_to_lineage, name_to_rank = tax_obj.load_taxonomy_map()

        random_colour_maker: RandomColour = RandomColour()

        def layout(node):
            nonlocal name_to_lineage, name_to_rank, random_colour_maker

            if node.is_leaf():
                if self.use_node_names:
                    faces.add_face_to_node(AttrFace("name", fsize=20), node, column=0, position="aligned")

                if self.at_rank is not None:
                    lineage = name_to_lineage[node.name]

                    for name in lineage:
                        rank = name_to_rank[name][1]

                        if rank == self.at_rank:
                            if name not in rank_colours:
                                rank_colours[name] = next(random_colour_maker)

                            node_style = NodeStyle()
                            node_style['bgcolor'] = rank_colours[name]
                            node.set_style(node_style)

                elif self.plot_phyla:
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

        phy_path = "phylotree.pdf" if not self.out_dir else os.path.join(self.out_dir, "phylotree.pdf")
        tree.render(
            phy_path,
            tree_style=ts,
            w=4000,
            units="px"
        )
        print("Phylogenetic tree written to %s!" % phy_path)

    """
    CLI for making phylogenetic tree
    ================================
    
    """
    def tree(self):
        conf: JSONConfig = self.get_conf()

        entry_points = (self.do_sketches, self.do_distances, self.do_trees)
        entry_stage = self.argmax(self.check_sketches(), self.check_distances(), self.check_trees())

        for stage in range(entry_stage, len(entry_points)):
            entry_points[stage]()

        tree_dir = self.finalise_tree()
        conf.set(phylogeny_path=tree_dir)
        conf.save()

        print("Phylogeny set to %s." % tree_dir)

    @staticmethod
    def argmax(*checks):
        if not any(checks):
            return 0

        return max(i for i, val in enumerate(checks + (True, )) if val is True)

    def finalise_tree(self):
        print("Finalising tree...")
        conf: JSONConfig = self.get_conf()

        tree_dir = os.path.join(self.config.phylogeny, 'tree')
        tree_name = "%s.nwk" % os.path.basename(self.config.phylogeny.rstrip(os.sep))
        tree_path = os.path.join(tree_dir, tree_name)

        tree_files = [
            name
            for name in ls(tree_dir, ext='.nwk')
            if not name.endswith(tree_name)
        ]

        def get_tree_data(file_url):
            with open(file_url, 'r') as f:
                tree_data = f.read().strip().rstrip(";")

            return tree_data

        if len(tree_files) == 1:
            nwk_data = get_tree_data(tree_files[0]).rstrip(";") + ";"

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

    """
    Make mashtree
    =============
    
    """
    def mashtree(self):
        self.do_mashtree()
        self.tree()

    def do_mashtree(self):
        print("Creating mashtree...")
        self.check_mashtree()
        conf: JSONConfig = self.get_conf()

        tree_dir = os.path.join(self.config.phylogeny, 'tree')
        temp_dir = os.path.join(self.config.phylogeny, 'tmp')
        n: int = self.get_n_processes()

        if not os.path.exists(tree_dir):
            os.mkdir(tree_dir)

        for group_name in conf.get_groups(self.group):
            k, s, sequences = self.get_group_or_die(group_name)
            tree_path = os.path.join(tree_dir, "%s.nwk" % group_name)

            self.make_mashtree(k, s, n, sequences, tree_path, temp_dir)

    def make_mashtree(self, k, s, n, sequences, tree_dir, temp_dir):
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

        except OSError:
            die("Could not write genome paths to text file!")

        # Make mashtree and insert into configuration file.
        print("Making mashtree...")
        cmd = "mashtree --numcpus %d --kmerlength %d --sketch-size %d --file-of-files %s" \
            % (n, k, s, _names_file)

        try:
            return_val = self.shell(cmd, cwd=temp_dir)

            with open(tree_dir, "w") as f:
                f.write(return_val)

            print("Phylogeny created in %s." % tree_dir)
        finally:
            try:
                delete_temp()
            except FileNotFoundError:
                pass
            except OSError:
                die("Could not delete temporary directory %s!" % temp_dir)

    """
    Sketch sequences
    ================
    
    """
    def sketch(self):
        self.do_sketches()

    def check_sketches(self):
        conf: JSONConfig = JSONConfig(self.config.conf)
        sketch_dir = os.path.join(self.config.phylogeny, 'sketch')
        sketch_name_fmt = "%s.k%d.s%d.%s"

        if not os.path.exists(sketch_dir):
            return False

        for group_name in conf.get_groups(self.group):
            k, s, _ = self.get_group_or_die(group_name)
            file_name = sketch_name_fmt % (group_name, k, s, "%s")

            dest = os.path.join(sketch_dir, file_name % ("sour" if self.use_sourmash else "msh"))
            if not os.path.exists(dest):
                return False

        return True

    def do_sketches(self):
        conf: JSONConfig = JSONConfig(self.config.conf)
        n: int = self.get_n_processes()

        sketch_dir = os.path.join(self.config.phylogeny, 'sketch')
        if not os.path.exists(sketch_dir):
            os.mkdir(sketch_dir)

        sketch_name_fmt = "%s.k%d.s%d.%s"

        for group_name in conf.get_groups(self.group):
            k, s, sequences = self.get_group_or_die(group_name)
            file_name = sketch_name_fmt % (group_name, k, s, "%s")

            if self.use_sourmash:
                sig_path = os.path.join(sketch_dir, file_name % "sour")
                if not os.path.exists(sig_path):        
                    os.mkdir(sig_path)

                self.sour_sketch(k=k, s=s, sequences=sequences, sig_dir=sig_path)

            else:
                self.check_mash()

                file_path = os.path.join(sketch_dir, file_name % "msh")
                self.mash_sketch(k=k, s=s, p=n, sequences=sequences, out_dir=file_path)

    def sour_sketch(self, k: int, s: int, sequences: List[str], sig_dir: str):
        from expam.tree.sourmash import make_signatures
        make_signatures(self.get_n_processes(), sequences, sig_dir, k, s)

    def mash_sketch(self, k, s, p, sequences, out_dir):
        cmd_fmt = "mash sketch -k %d -p %d -s %d -o %s %s"
        cmd = cmd_fmt % (k, p, s, out_dir, ' '.join(sequences))
        self.shell(cmd)

    """
    Compute pairwise distances between sketched sequences
    =====================================================
    
    """
    def distance(self):
        self.do_distances()

    def check_distances(self):
        conf: JSONConfig = self.get_conf()

        sketch_dir = os.path.join(self.config.phylogeny, 'sketch')
        dist_dir = os.path.join(self.config.phylogeny, 'distance')

        if not os.path.exists(dist_dir):
            return False

        dist_name_fmt = "%s.k%d.s%d.tab"

        for group_name in conf.get_groups(self.group):
            k, s, _ = self.get_group_or_raise(group_name)
            file_name = dist_name_fmt % (group_name, k, s)

            dest = os.path.join(sketch_dir, file_name)
            if not os.path.exists(dest):
                return False

        return True

    def do_distances(self):
        print("Calculating pairwise distances...")
        conf: JSONConfig = self.get_conf()

        sketch_dir = os.path.join(self.config.phylogeny, 'sketch')
        dist_dir = os.path.join(self.config.phylogeny, 'distance')

        if not os.path.exists(dist_dir):
            os.mkdir(dist_dir)

        dist_name_fmt = "%s.k%d.s%d.%s"

        for group_name in conf.get_groups(self.group):
            k, s, _ = self.get_group_or_die(group_name)
            sketch_name = dist_name_fmt % (group_name, k, s, '%s')

            matrix_name = dist_name_fmt % (group_name, k, s, 'tab')
            matrix_file = os.path.join(dist_dir, matrix_name)

            if self.use_sourmash:
                sketch_dir = os.path.join(sketch_dir, sketch_name % 'sour')
                self.sour_dist(sig_dir=sketch_dir, matrix_dir=matrix_file)

            else:
                self.check_mash()

                sketch_file = os.path.join(sketch_dir, sketch_name % 'msh')
                self.mash_dist(sketch_file, matrix_file)

    def mash_dist(self, sketch_dir, matrix_dir):
        self.check_mash()

        cmd_fmt = "mash dist -p %d -t %s %s"
        cmd = cmd_fmt % (
            self.get_n_processes(),
            sketch_dir,
            sketch_dir
        )

        _unfmt_dst = self.shell(cmd)
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

    def sour_dist(self, sig_dir, matrix_dir):
        from expam.tree.sourmash import make_distances
        make_distances(sig_dir, self.get_n_processes(), matrix_dir)

    """
    Execute neighbour-joining on distance matrix
    ============================================
    
    """
    def nj(self):
        self.do_trees()

    def check_trees(self):
        conf: JSONConfig = self.get_conf()

        tree_dir = os.path.join(self.config.phylogeny, 'tree')
        if not os.path.exists(tree_dir):
            return False

        for group_name in conf.get_groups(self.group):
            tree_name = '%s.nwk' % group_name
            tree_path = os.path.join(tree_dir, tree_name)

            if not os.path.exists(tree_path):
                return False

        return True

    def do_trees(self):
        print("Running neighbour-joining...")
        conf: JSONConfig = self.get_conf()

        dist_dir = os.path.join(self.config.phylogeny, 'distance')
        tree_dir = os.path.join(self.config.phylogeny, 'tree')
        dist_fmt = "%s.k%d.s%d.tab"

        if not os.path.exists(tree_dir):
            os.mkdir(tree_dir)

        for group_name in conf.get_groups(self.group):
            k, s, _ = self.get_group_or_die(group_name)

            matrix_name = dist_fmt % (group_name, k, s)
            matrix_path = os.path.join(dist_dir, matrix_name)

            tree_name = "%s.nwk" % group_name
            tree_path = os.path.join(tree_dir, tree_name)

            if self.use_quicktree:
                self.quicktree(matrix_dir=matrix_path, tree_dir=tree_path)
            else:
                self.rapidnj(matrix_dir=matrix_path, tree_dir=tree_path)

    def rapidnj(self, matrix_dir, tree_dir):
        self.check_rapidnj()

        cmd_fmt = "rapidnj %s -i pd -o t -c %d"
        cmd = cmd_fmt % (matrix_dir, self.get_n_processes())

        unformatted_tree = self.shell(cmd)
        tree = self._format_tree_string(unformatted_tree)

        with open(tree_dir, 'w') as f:
            f.write(tree)

        print("RapidNJ wrote tree to %s." % tree_dir)

    def quicktree(self, matrix_dir, tree_dir):
        self.check_quicktree()

        cmd_fmt = "quicktree -in m -out t %s"
        cmd = cmd_fmt % matrix_dir

        unformatted_tree = self.shell(cmd)
        tree = self._format_tree_string(unformatted_tree)

        with open(tree_dir, 'w') as f:
            f.write(tree)

        print("QuickTree wrote tree to %s." % tree_dir)

    @staticmethod
    def _format_tree_string(nwk):
        seq_types = [".fna", ".faa"]
        comp_types = [".tar.gz", ".gz"]

        nwk = nwk.replace("'", "")

        for seq_type in seq_types:
            for comp_type in comp_types:
                ext = seq_type + comp_type + ":"
                nwk = nwk.replace(ext, ":")

        return nwk
