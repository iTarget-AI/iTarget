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


cd ./_ForFeatures/esm2
conda activate iTarget
cp ./bashes/template_esm2_t36_3B_UR50D.sh ./bashes/template_esm2_t36_3B_UR50D.runtmp.sh
sed -i "s|--datatype template|--datatype ${PROTROOT}|g" ./bashes/template_esm2_t36_3B_UR50D.runtmp.sh
cd ./bashes
bash template_esm2_t36_3B_UR50D.runtmp.sh
cd ..
# sed -i "s|--datatype ${PROTROOT}|--datatype template|g" ./bashes/template_esm2_t36_3B_UR50D.sh
cp ./data/esm2_t36_3B_UR50D/${PROTROOT}/${PROTROOT}_all-data-merge-prot.csv ../../data/original_data
cd ../..


cd ./_ForFeatures/xmol
conda activate iTarget
cp ./bashes/template_xmol.sh ./bashes/template_xmol.runtmp.sh
sed -i "s|self|${DRUGROOT}|g" ./bashes/template_xmol.runtmp.sh
cd ./bashes
bash template_xmol.runtmp.sh
cd ..
# sed -i "s|${DRUGROOT}|self|g" ./bashes/template_xmol.sh
cp ./FT_to_embedding/data/for_output/${DRUGROOT}_all-data-merge-drug.csv ../../data/original_data
cd ../..
