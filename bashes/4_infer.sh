#!/usr/bin/bash

python -u ../scripts/main.py --kfold_num 1 --task predict --gpu 0 --batch_size 1 --source infer --pretrained_model_path pretrained/iTarget_final_model --uncertainty_iteration_size 100
