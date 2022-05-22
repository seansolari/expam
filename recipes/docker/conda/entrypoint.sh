#!/bin/bash --login
# The --login ensures the bash configuration is loaded,
# enabling Conda.

# Temporarily disable strict mode and activate conda:
set +euo pipefail
conda activate expam

# Re-enable strict mode:
set -euo pipefail

expam --help
