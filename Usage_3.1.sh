#!/usr/bin/bash
set -euo pipefail
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate iTarget


# export PROTROOT=self
# export DRUGROOT=self
export PROTROOT=bindingdb
export DRUGROOT=bindingdb


cd ./bashes
# for compounds, by default, {--scale_method}='standard', {--disttype}='uniprot+fullchembl', {--source}='example' or user-defined
cp ./1_trans_drug.sh ./1_trans_drug.runtmp.sh
sed -i "s|--source example|--source ${DRUGROOT}|g" ./1_trans_drug.runtmp.sh
bash 1_trans_drug.runtmp.sh

# for proteins, by default, {--scale_method}='standard', {--disttype}='uniprot+fullchembl', {--source}='example' or user-defined
cp ./1_trans_prot.sh ./1_trans_prot.runtmp.sh
sed -i "s|--source example|--source ${PROTROOT}|g" ./1_trans_prot.runtmp.sh
bash 1_trans_prot.runtmp.sh


# navigate back to the project root directory
cd .. 