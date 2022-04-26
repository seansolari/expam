import os

from expam.classify import PHY_CLASSIFIED_FILE, PHY_RAW, PHY_RESULTS, PHY_SPLIT_FILE, TAX_CLASSIFIED_FILE, TAX_RAW, TAX_RESULTS, TAX_SPLIT_FILE, TEMP_RESULTS, ResultsPathConfig
from expam.utils import die


def make_results_config(out_path: str) -> ResultsPathConfig:
    output_file_locations = {
        'base': out_path,
        'phy': os.path.join(out_path, PHY_RESULTS),
        'tax': os.path.join(out_path, TAX_RESULTS),
        'temp': os.path.join(out_path, TEMP_RESULTS),
        'phy_raw': os.path.join(out_path, PHY_RAW),
        'tax_raw': os.path.join(out_path, TAX_RAW),
        'phy_classified': os.path.join(out_path, PHY_CLASSIFIED_FILE),
        'phy_split': os.path.join(out_path, PHY_SPLIT_FILE),
        'tax_classified': os.path.join(out_path, TAX_CLASSIFIED_FILE),
        'tax_split': os.path.join(out_path, TAX_SPLIT_FILE)
    }
    return ResultsPathConfig(**output_file_locations)

def load_results_config(out_path: str, create: bool = False) -> ResultsPathConfig:
    proposed_config: ResultsPathConfig = make_results_config(out_path)

    # Make base results path.
    if create:
        if not os.path.exists(out_path):
            try:
                os.mkdir(out_path)
            except OSError:
                print("Failed to make results path %s." % out_path)

        create_results(proposed_config)
    elif not validate_results_configuration(proposed_config, check_taxonomy=False):
        die("Results path does not exist!")
        
    return proposed_config

def create_results(config: ResultsPathConfig):
    for path_field in ('phy', 'tax', 'phy_raw', 'tax_raw', 'temp'):
        path = getattr(config, path_field)

        if not os.path.exists(path):
            os.mkdir(path)

def validate_results_configuration(config: ResultsPathConfig, check_taxonomy: bool = True):
    phy_files = (config.phy, config.phy_classified, config.phy_split)
    tax_files = (config.tax, config.tax_classified, config.tax_split)

    for phy_file in phy_files:
        if not os.path.exists(phy_file):
            return False

    if check_taxonomy:
        for tax_file in tax_files:
            if not os.path.exists(tax_file):
                return False
    
    return True

def validate_classification_results(results_dir: str):
    if not os.path.exists(results_dir):
        die("Could not find results %s!" % results_dir)

    results_config: ResultsPathConfig = make_results_config(results_dir)

    if not (os.path.exists(results_config.phy_classified) or os.path.exists(results_config.phy_split)):
        raise Exception("Path does not look like expam results folder!")

