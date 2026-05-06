#!/usr/bin/bash
set -euo pipefail

# This repository tracks large assets as URL-backed DVC imports rather than a
# configured DVC remote. Use `dvc update` to fetch them on a new machine.

dvc update _conda_envs/iTarget.tar.gz.dvc
dvc update _ForFeatures/xmol/FT_to_embedding/python.tar.gz.dvc
dvc update _ForFeatures/xmol/FT_to_embedding/data/model/step_400000/step_400000_20200326221400.tar.dvc
DVC_HTTP_TIMEOUT=3600 dvc update _ForFeatures/esm2/pretrained_esm2_models/esm2_t36_3B_UR50D.pt.dvc
dvc update data/processed_data/esm2_fitted.mp.dvc
dvc update data/_original_data/chembl_cpi/ChEMBL_data.csv.dvc
