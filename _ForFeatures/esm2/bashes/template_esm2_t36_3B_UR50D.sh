#!/usr/bin/bash
set -euo pipefail
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate iTarget

# esm2type = 'esm2_t33_650M_UR50D'
# repr_layer = 33
# hid_dim = 1280
# esm2type = 'esm2_t36_3B_UR50D'
# repr_layer = 36
# hid_dim = 2560
# esm2type = 'esm2_t48_15B_UR50D'
# repr_layer = 48
# hid_dim = 5120

echo "ESM2 is Running......May Take Hours......Please Wait......"
python -u ../main.py --start_index head --esm2type esm2_t36_3B_UR50D --datatype template --repr_layer 36 --hid_dim 2560
