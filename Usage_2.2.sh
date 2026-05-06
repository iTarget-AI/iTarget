#!/usr/bin/bash
set -euo pipefail
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate iTarget


# export PROTROOT=uniprot
# export DRUGROOT=fullchembl
export PROTROOT=self
export DRUGROOT=self
# export PROTROOT=bindingdb
# export DRUGROOT=bindingdb

export SCALE_COURCE="${PROTROOT}+${DRUGROOT}"


# calculate feature distance
cp ./bashes/0_feadist.sh ./bashes/0_feadist.runtmp.sh
sed -i "s|--scale_source uniprot+fullchembl|--scale_source ${SCALE_COURCE}|g" ./bashes/0_feadist.runtmp.sh
cd ./bashes
bash 0_feadist.runtmp.sh	# by default, {--scale_method}='standard', {--scale_source}='uniprot+fullchembl'


# copy calculated configs to work path
mkdir -p ../feamap/config/trans_from_${SCALE_COURCE}
cp ../data/processed_data/drug_fea/scale/standard/*.cfg ../feamap/config/trans_from_${SCALE_COURCE}/
cp ../data/processed_data/protein_fea/scale/standard/*.cfg ../feamap/config/trans_from_${SCALE_COURCE}/

# navigate back to the project root directory
cd .. 