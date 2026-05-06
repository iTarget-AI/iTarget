#!/usr/bin/bash
set -euo pipefail
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate iTarget


# export BENCHMARK=human
# export BENCHMARK=biosnap
# export BENCHMARK=chembl

# optional, or you can directly prepare files following the examples in `./data/processed_data/split_cvdata/`.
# This step is not required for bindingdb benchmark, which has been done in step 1.1
cd ./bashes
cp ./2_split_cvdata.sh ./2_split_cvdata.runtmp.sh
sed -i "s|--source example|--source ${BENCHMARK}|g" ./2_split_cvdata.runtmp.sh
bash 2_split_cvdata.runtmp.sh

# navigate back to the project root directory
cd .. 