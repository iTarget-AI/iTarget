#!/usr/bin/bash
set -euo pipefail
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate iTarget


export ROOT=bindingdb
# export ROOT=human
# export ROOT=biosnap
# export ROOT=chembl

cd ./bashes
cp ./3_train_cv.sh ./3_train_cv.runtmp.sh
sed -i "s|--source example|--source ${ROOT}|g" ./3_train_cv.runtmp.sh
bash 3_train_cv.runtmp.sh # by defalut, {--kfold_num}=5, {--task}='cv', {--n_epochs}=128, {--gpu}=0, {--batch_size}=512, {--lr}=5e-4, {--monitor}='auc_val', {--source}='example', {--uncertainty_iteration_size}=100


# navigate back to the project root directory
cd .. 