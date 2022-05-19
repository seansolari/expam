import os
from typing import List, Set
from expam.cli.main import CommandGroup, ExpamOptions, clear_logs
from expam.database import FileLocationConfig
from expam.database.build import main as expam
from expam.database.config import ExpamDatabaseDoesNotExistError, JSONConfig, create_database, make_database_config, validate_database_file_configuration
from expam.logger import Timer
from expam.utils import die, ls, make_path_absolute


class BuildCommand(CommandGroup):
    commands: Set[str] = {
        'quickrun', 'default_db', 'create',
        'build', 'print', 'add', 'remove', 'set'
    }

    def __init__(
        self, config: FileLocationConfig,
        k: int, n: int, s: int, phylogeny_path: str, pile_size: int,
        files: List[str], group: str, first_n: int
    ) -> None:
        super().__init__()
        self.config: FileLocationConfig = config

        self.k = k
        self.n = n
        self.s = s
        self.phylogeny_path = phylogeny_path
        self.pile_size = pile_size

        self.files = files
        self.group = group
        self.first_n = first_n

    @classmethod
    def take_args(cls: CommandGroup, args: ExpamOptions) -> dict:
        k, n, s, pile_size, first_n = cls.parse_ints(args.k, args.n, args.s, args.pile, args.first_n)

        return {
            'config': make_database_config(args.db_name),
            'k': k,
            'n': n,
            's': s,
            'phylogeny_path': args.phylogeny,
            'pile_size': pile_size,
            'files': args.directory,
            'group': None if args.groups is None else args.groups[0][0],
            'first_n': first_n
        }

    def validate_database(self):
        try:
            validate_database_file_configuration(self.config)
        except ExpamDatabaseDoesNotExistError:
            die("%s is not a database." % self.config.database)

    """
    Quickrun command
    ================
    
    """
    def quickrun(self):
        self.create()
        self.set()
        self.add()
        self.build()

    """
    Default database command
    ========================
    
    """
    def default_db(self):
        self.validate_database()

        print("export EXPAM_DB_DIR=%s" % self.config.database)

    """
    Create database command
    =======================
    
    """
    def create(self):
        create_database(self.config)

    """
    Build database command
    ======================
    
    """
    def build(self):
        self.validate_database()
        conf = JSONConfig(self.config.conf)

        # Check if a phylogeny has been provided.
        if not self.check_phylogeny(conf):
            die("No phylogeny set for this database.")

        # Read configuration file.
        k, n, phylogeny, genome_paths, pile = conf.get_build_params()

        clear_logs(self.config.logs)

        # Set building database.
        try:
            with Timer() as t:
                expam(
                    db_path=self.config.base,
                    genome_paths=genome_paths,
                    phylogeny_path=phylogeny,
                    k=k,
                    n=n - 1,  # Account for main process.
                    pile_size=pile,
                )

            print("expam: " + str(t))

        finally:
            print(f'\nPID - {os.getpid()} dying...')

    @staticmethod
    def check_phylogeny(conf):
        if conf['phylogeny_path'] is None:
            return False
        else:
            return True

    """
    Print database command
    ======================
    
    """
    def print(self):
        self.validate_database()

        conf: JSONConfig = JSONConfig(self.config.conf)
        print(conf)

    """
    Add sequences to database command
    =================================
    
    """
    def add(self):
        self.validate_database()

        for directory in self.files:
            self.add_sequences(directory)

    def add_sequences(self, path):
        self._modify_config(path, add=True)

    def _modify_config(self, directory, add=True):
        base_url = make_path_absolute(
            dir=directory,
            base_dir=os.getcwd()
        )

        # Get list of sequences in base_url.
        seq_names = ls(base_url)

        if (self.first_n is None) or (self.first_n > len(seq_names)):
            self.first_n = len(seq_names)

        seq_dirs = [
            os.path.join(base_url, seq_names[i])
            if seq_names[i] != base_url else seq_names[i]
            for i in range(self.first_n)
        ]

        # Append directory to database configuration.
        conf: JSONConfig = JSONConfig(self.config.conf)

        if add:
            conf.add_sequence(seq_dirs, group=self.group)
            print("Added %d files from %s." % (len(seq_dirs), base_url))
        else:
            conf.remove_sequence(seq_dirs, group=self.group)
            print("Removed %d files from %s." % (len(seq_dirs), base_url))

        conf.save()

    """
    Remove sequences to database command
    ====================================
    
    """
    def remove(self):
        self.validate_database()

        for directory in self.files:
            self.add_sequences(directory)

    def remove_sequences(self, path):
        self._modify_config(path, add=False)

    
    """
    Set database settings
    =====================
    
    """
    def set(self):
        self.validate_database()

        conf = JSONConfig(url=self.config.conf)

        if self.group is None:
            conf.set(k=self.k, n=self.n, s=self.s, phylogeny_path=self.phylogeny_path, pile=self.pile_size)
        else:
            conf.group_set(group_name=self.group, k=self.k, s=self.s)

        conf.save()
