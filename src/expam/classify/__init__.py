from collections import namedtuple

import os


ENTREZ_EMAIL = "sean.solari@monash.edu"
ENTREZ_TOOL = "expam"
ENTREZ_API = "8e5edfa1413576bca4f48b4a5520f6295308"
_STEP = 500  # NCBI Database searches at a time.

PHY_RESULTS = "phy"
TAX_RESULTS = "tax"
RAW_RESULTS = "raw"
TEMP_RESULTS = "temp"

CLASSIFIED_NAME = "classified.tsv"
SPLIT_NAME = "split.tsv"

PHY_RAW = os.path.join(PHY_RESULTS, RAW_RESULTS)
TAX_RAW = os.path.join(TAX_RESULTS, RAW_RESULTS)
PHY_CLASSIFIED_FILE = os.path.join(PHY_RESULTS, CLASSIFIED_NAME)
PHY_SPLIT_FILE = os.path.join(PHY_RESULTS, SPLIT_NAME)
TAX_CLASSIFIED_FILE = os.path.join(TAX_RESULTS, CLASSIFIED_NAME)
TAX_SPLIT_FILE = os.path.join(TAX_RESULTS, SPLIT_NAME)


ResultsPathConfig = namedtuple(
    'ResultsPathConfig',
    [
        'base', 'phy', 'tax',
        'temp',
        'phy_raw', 'tax_raw',
        'phy_classified', 'phy_split',
        'tax_classified', 'tax_split'
    ]
)
