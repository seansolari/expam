import json
import os

import numpy as np
from expam.database import ACCESSION_ID_RELATIVE_PATH, CONF_RELATIVE_PATH, DATABASE_FILE_RELATIVE_PATH, DATABASE_RELATIVE_PATH, LCA_MATRIX_RELATIVE_PATH, LOG_RELATIVE_PATH, PHYLOGENY_RELATIVE_PATH, TAXID_LINEAGE_MAP_RELATIVE_PATH, TAXON_RANK_MAP_RELATIVE_PATH, FileLocationConfig
from expam.utils import die


class JSONConfig:
    _parse_data_type = {
        "uint8": np.uint8,
        "uint16": np.uint16,
        "uint32": np.uint32,
        "uint64": np.uint64,
        "int8": np.int8,
        "int16": np.int16,
        "int32": np.int32,
        "int64": np.int64,
    }

    def __init__(self, url=None, name: str = ""):
        self.url = url

        # Default values.
        self.params = {
            "name": name,
            "groups": {
                'default': {
                    'sequence_files': [],
                    'phylogeny_path': [],
                    'k': None,
                    's': None,
                },
            },
            "phylogeny_path": None,
            "k": None,
            "n": None,
            "s": None,
            "pile": None,
            "keys_shape": None,
            "keys_data_type": None,
            "values_shape": None,
            "values_data_type": None,
        }

        if url is not None:
            with open(url, 'r') as f:
                _conf = json.load(f)

                if not isinstance(_conf, dict):
                    die('Invalid configuration file!')

                for key in _conf:
                    if key == 'groups':
                        self.params['groups'].update(_conf['groups'])

                    else:
                        self.params[key] = _conf[key]

    def __getitem__(self, key):
        return self.params[key]

    def __str__(self):
        name_line = "<<< expam configuration file: %s >>>\n" % self.params['name']
        data_line = "phylogeny\t-->\t%s\nk\t\t-->\t%s\nn\t\t-->\t%s" % (self.params['phylogeny_path'],
                                                                    self.params['k'], self.params['n'])
        par_line = "sketch\t\t-->\t%s\npile\t\t-->\t%s" % (self.params['s'], self.params['pile'])

        group_fmt = "\n----------------\ngroup name: %s\n\tk\t\t-->\t%s\n\tsketch\t\t-->\t%s\n\tsequences\t-->\t%d\n"
        group_line = "\n".join([
            group_fmt % (name, group['k'], group['s'], len(group['sequence_files']))
            for name, group in self.params['groups'].items()
        ])

        return "\n".join([name_line, data_line, par_line, group_line])

    def get_paths(self):
        return [
            genome_path
            for group in self.params["groups"].values()
            for genome_path in group['sequence_files']
        ]

    def new_group(self, name, files=None, phylogeny_path=None, k=None, s=None):
        if name in self.params['groups']:
            raise ValueError("Group %s already exists!" % name)

        self.params['groups'][name] = {
            'sequence_files': [] if files is None else files,
            'phylogeny_path': phylogeny_path,
            'k': k,
            's': s,
        }

    def set(self, **kwargs):
        for arg in kwargs:
            if arg in self.params and arg != 'groups':
                if kwargs[arg] is not None:
                    self.params[arg] = kwargs[arg]

            else:
                raise ValueError("Unknown configuration category %s!" % arg)

    def group_set(self, group_name, **kwargs):
        for arg in kwargs:
            if arg in self.params['groups'][group_name] and arg != 'sequence_files' and kwargs[arg] is not None:
                self.params['groups'][group_name][arg] = kwargs[arg]

    def group_get(self, group_name):
        if group_name not in self.params['groups']:
            raise IndexError("Invalid group %s!" % group_name)

        group = self.params['groups'][group_name]
        attrs = {attr: group[attr] for attr in ('k', 's', 'sequence_files')}

        for req_param in ('k', 's'):
            if group[req_param] is None:
                if self.params[req_param]:
                    attrs[req_param] = int(self.params[req_param])

            else:
                attrs[req_param] = int(attrs[req_param])

        return attrs['k'], attrs['s'], attrs['sequence_files']

    def add_sequence(self, file_name, group=None):
        if group is None:
            group = 'default'

        if group not in self.params['groups']:
            self.new_group(group)

        if not isinstance(file_name, list):
            file_name = {file_name}
        else:
            file_name = set(file_name)

        # Check these sequences aren't already added.
        for current_group in self.params['groups']:
            if current_group != group:
                current_seqs = set(self.params['groups'][current_group]['sequence_files'])

                if file_name & current_seqs:
                    die("Some of these sequences are already in the database!")

        self.params['groups'][group]['sequence_files'].extend(file_name)

    def remove_sequence(self, file_name, group=None):
        if group is None:
            group = 'default'

        if isinstance(file_name, list):
            file_name = set(file_name)
        else:
            file_name = {file_name}

        sequences = set(self.params['groups'][group]['sequence_files'])
        sequences -= file_name

        self.params['groups'][group]['sequence_files'] = list(sequences)

    def remove_group(self, group=None):
        if group is None:
            self.params['groups']['default'] = {
                    'sequence_files': [],
                    'phylogeny_path': [],
                    'k': None,
                    's': None,
                }
        elif group in self.params['groups']:
            self.params['groups'].pop(group)
        else:
            raise Exception("No such group %s!" % group)

    def _dump(self):
        return json.dumps(self.params)

    def groups(self):
        return list(self.params['groups'].keys())

    def save(self, url=None):
        if url is None:
            if self.url is None:
                raise Exception("Save url required!")

            url = self.url

        with open(url, 'w') as f:
            json.dump(self.params, f)

    def get_build_params(self):
        """Returns k, n, phylogeny_path, genome_paths, pile_size"""
        genome_paths, phylogeny_path, k, n, pile_size = self.get_paths(), self["phylogeny_path"], int(self["k"]), int(self["n"]), self["pile"]
        return k, n, phylogeny_path, genome_paths, pile_size

    def get_groups(self, input_group: str = None):
        groups = self.groups() if input_group is None else [input_group]

        # Check group has more than zero sequences.
        for group in groups:
            _, _, sequences = self.group_get(group)

            if len(sequences) > 0:
                yield group

    def get_n_processes(self):
        return self['n']


class ExpamDatabaseDoesNotExistError(Exception):
    pass

class ExpamDatabaseExistsError(Exception):
    pass


def make_database_config(db_path: str) -> FileLocationConfig:
    database_file_locations = {
        'base': db_path,
        'database': os.path.join(db_path, DATABASE_RELATIVE_PATH),
        'phylogeny': os.path.join(db_path, PHYLOGENY_RELATIVE_PATH),
        'logs': os.path.join(db_path, LOG_RELATIVE_PATH),
        'conf': os.path.join(db_path, CONF_RELATIVE_PATH),
        'accession_id': os.path.join(db_path, ACCESSION_ID_RELATIVE_PATH),
        'taxid_lineage': os.path.join(db_path, TAXID_LINEAGE_MAP_RELATIVE_PATH),
        'taxon_rank': os.path.join(db_path, TAXON_RANK_MAP_RELATIVE_PATH),
        'lca_matrix': os.path.join(db_path, LCA_MATRIX_RELATIVE_PATH),
        'database_file': os.path.join(db_path, DATABASE_FILE_RELATIVE_PATH),
    }
    return FileLocationConfig(**database_file_locations)


def load_database_config(db_path: str) -> FileLocationConfig:
    proposed_config: FileLocationConfig = make_database_config(db_path)

    validate_database_file_configuration(proposed_config)
    return proposed_config


def create_database(config: FileLocationConfig) -> None:
    for field_to_check in ('base', 'database', 'phylogeny', 'logs'):
        path: str = getattr(config, field_to_check)

        if not os.path.exists(path):
            os.mkdir(path)
        else:
            raise ExpamDatabaseExistsError("Database %s already exists!" % config.database)

    # Create new configuration file.
    conf: JSONConfig = JSONConfig()
    conf.save(url=config.conf)


def validate_database_file_configuration(proposed_config: FileLocationConfig) -> bool:
    for field_to_check in ('base', 'database', 'phylogeny', 'conf'):
        path_to_check = getattr(proposed_config, field_to_check)

        if not os.path.exists(path_to_check):
            print(path_to_check)
            raise ExpamDatabaseDoesNotExistError("Database does not exist at %s" % path_to_check)


def validate_taxonomy_files(config: FileLocationConfig) -> bool:
    for field_to_check in ('accession_id', 'taxid_lineage', 'taxon_rank'):
        if not os.path.exists(getattr(config, field_to_check)):
            return False
    else:
        return True
