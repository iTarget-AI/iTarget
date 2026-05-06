#!/usr/bin/bash
set -euo pipefail
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate iTarget


cd ./data/predict_data/infer
cp prepare_predict_files.py prepare_predict_files.runtmp.py
python prepare_predict_files.runtmp.py
cd ..
mv infer_prots.csv ../../_ForFeatures/esm2/data
mv infer_drugs.csv ../../_ForFeatures/xmol/FT_to_embedding/data/for_output
cd ../..


cd ./_ForFeatures/esm2
cp ./bashes/template_esm2_t36_3B_UR50D.sh ./bashes/template_esm2_t36_3B_UR50D.runtmp.sh
sed -i "s|--datatype template|--datatype infer|g" ./bashes/template_esm2_t36_3B_UR50D.runtmp.sh
cd ./bashes
bash template_esm2_t36_3B_UR50D.runtmp.sh
cd ..
cp ./data/esm2_t36_3B_UR50D/infer/infer_all-data-merge-prot.csv ../../data/predict_data


cd ../xmol
cp ./bashes/template_xmol.sh ./bashes/template_xmol.runtmp.sh
sed -i "s|self|infer|g" ./bashes/template_xmol.runtmp.sh
cd ./bashes
bash template_xmol.runtmp.sh
cd ..
cp ./FT_to_embedding/data/for_output/infer_all-data-merge-drug.csv ../../data/predict_data


cd ..
