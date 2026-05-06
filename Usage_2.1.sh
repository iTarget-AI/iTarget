#!/usr/bin/bash

# export PROTROOT=uniprot
# export DRUGROOT=fullchembl
export PROTROOT=self
export DRUGROOT=self
# export PROTROOT=bindingdb
# export DRUGROOT=bindingdb

export SCALE_COURCE="${PROTROOT}+${DRUGROOT}"

# cp ./_ForFeatures/esm2/data/{--esm2type}/{--datatype}/{--datatype}_all-data-merge-prot.csv ./data/original_data/scale/
# for proteins' template, by default, {--esm2type}='esm2_t36_3B_UR50D', {--datatype}='uniprot'
cp ./_ForFeatures/esm2/data/esm2_t36_3B_UR50D/${PROTROOT}/${PROTROOT}_all-data-merge-prot.csv ./data/original_data/scale/

# mv ./_ForFeatures/xmol/FT_to_embedding/data/for_output/{--datatype}_all-data-merge-drug.csv ./data/original_data/scale/
# for compounds' template, by default, {--datatype}='fullchembl'
cp ./_ForFeatures/xmol/FT_to_embedding/data/for_output/${DRUGROOT}_all-data-merge-drug.csv ./data/original_data/scale/


# The moved feature files in `./data/original_data/scale/`  should be renamed using same {--scale_source} for {--datatype} according to the corresponding settings in downstream file `0_feadist.sh`.
# Here, we use 'uniprot+fullchembl' as an example, and then result in `uniprot+fullchembl_all-data-merge-prot.csv` and `uniprot+fullchembl_all-data-merge-drug.csv` two files in `./data/original_data/scale/`.
mv ./data/original_data/scale/${PROTROOT}_all-data-merge-prot.csv ./data/original_data/scale/${SCALE_COURCE}_all-data-merge-prot.csv
mv ./data/original_data/scale/${DRUGROOT}_all-data-merge-drug.csv ./data/original_data/scale/${SCALE_COURCE}_all-data-merge-drug.csv
