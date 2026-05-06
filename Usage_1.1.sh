#!/usr/bin/bash
set -euo pipefail
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate iTarget

python scripts/_data_preprocess.py
# the produced '{type}_drugs.csv' and '{type}_prots.csv' files could be used in step 1.2
