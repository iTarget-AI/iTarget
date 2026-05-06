#!/usr/bin/bash
set -euo pipefail
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate iTarget


cd ./bashes
cp ./4_infer.sh ./4_infer.runtmp.sh
sed -i "s|--source infer|--source infer|g" ./4_infer.runtmp.sh
bash 4_infer.runtmp.sh # by defalut, {--kfold_num}=1, {--task}='predict', {--gpu}=0, {--batch_size}=1, {--source}='infer', {--pretrained_model_path}='pretrained/iTarget_final_model', {--uncertainty_iteration_size}=100

# navigate back to the project root directory
cd .. 