from collections import namedtuple
import logging
import os
import numpy as np

KEYS_DATA_TYPE = np.uint64
KEYS_DATA_TYPE_SIZE = np.dtype(KEYS_DATA_TYPE).itemsize
KEYS_DTYPE_STR = "uint64"

MAP_DATA_TYPE = np.uint16
MAP_DATA_TYPE_SIZE = np.dtype(MAP_DATA_TYPE).itemsize
MAP_DTYPE_STR = "uint16"

DataTypeConfig = namedtuple('DataTypeConfig', ['keys_dtype', 'keys_dtype_size', 'keys_dtype_str', 'map_dtype', 'map_dtype_size', 'map_dtype_str'])
expam_dtypes = DataTypeConfig(KEYS_DATA_TYPE, KEYS_DATA_TYPE_SIZE, KEYS_DTYPE_STR, MAP_DATA_TYPE, MAP_DATA_TYPE_SIZE, MAP_DTYPE_STR)

DATABASE_RELATIVE_PATH = "database"
PHYLOGENY_RELATIVE_PATH = "phylogeny"
LOG_RELATIVE_PATH = "logs"
CONF_RELATIVE_PATH = "conf.json"

# There are five taxonomy data files:
#   - ACCN_ID       ... Maps genome sequences to taxonomy IDs.
ACCN_ID = os.path.join(PHYLOGENY_RELATIVE_PATH, "accession_ids.csv")
#   - CLADE_TAX     ... Maps phylogenetic clades to taxonomy IDs - depends on ACCN_ID.
CLADE_TAX = os.path.join(PHYLOGENY_RELATIVE_PATH, "clade_taxids.csv")
#   - TAX_POOL      ... Lists taxonomy IDs such that parent taxids always occur before children.
TAX_POOL = os.path.join(PHYLOGENY_RELATIVE_PATH, "taxid_pool.txt")
#   - TAX_CHILDS    ... Maps taxids to immediate child node taxids.
TAX_CHILDS = os.path.join(PHYLOGENY_RELATIVE_PATH, "child_taxids.csv")
#   - CLADE_TABLE   ... Lists taxonomy IDs for clades, with taxonomic data
CLADE_TABLE = os.path.join(PHYLOGENY_RELATIVE_PATH, "clade_taxonomy.csv")

LCA_MATRIX_RELATIVE_PATH = os.path.join(PHYLOGENY_RELATIVE_PATH, "lca_matrix.npy")
DATABASE_FILE_RELATIVE_PATH = os.path.join(DATABASE_RELATIVE_PATH, "expam_db.h5")

FileLocationConfig = namedtuple(
    'FileLocationConfig',
    [
        'base', 'database', 'phylogeny', 'logs', 'conf',
        'accn_id', 'clade_tax', 'tax_pool', 'tax_childs', 'clade_table',
        'lca_matrix',
        'database_file'
    ]
)

CHUNK_SIZE = 32
UNION_RATIO = 4 / 5

TIMEOUT = 1e-4

NULL_VALUE = 0  # Return value for foreign kmer.
EXPAM_TEMP_EXT = "expam_temp"

LOG_FORMAT = logging.Formatter(
    fmt='%(asctime)s... %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p'
)
