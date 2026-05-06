#!/usr/bin/bash

python -u ../scripts/0_feadist_fit.py --overwrite_cache --scale_method standard --drug_featype xmol --prot_featype esm2 --scale_source uniprot+fullchembl
