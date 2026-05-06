# iTarget



![iTarget.png](iTarget.png)



## Installation

##### 1. Platform Requirements

iTarget is designed to run on Linux, requiring the prior installation of `Anaconda`, `Git`, and `DVC`.

```
https://www.anaconda.com/docs/getting-started/anaconda/install/overview
https://git-scm.com/install
https://doc.dvc.org/install
```

##### 2. Download the Source Code and Required Resources

This repository contains the full source code of the project. Because Git does not track large assets efficiently, pretrained model weights, large datasets, and packaged Conda environments are not fully stored in the **Git** repository.

Clone the repository with Git:

```bash
git clone git@github.com:iTarget-AI/iTarget.git
cd ./iTarget
```

This repository stores large assets as URL-backed DVC imports. On a new machine,
use **DVC** to fetch those imported files from their recorded download links:

```bash
bash bashes/pull_external_assets.sh
```

**Alternatively**, **if preferred or if the automatic download fails**, you may download the required files manually and place them in the corresponding directories listed below.

| File Name | Description | File Size | Target Directory | Download Link |
| --- | --- | --- | --- | --- |
| `iTarget.tar.gz` | Packaged runtime environment for iTarget | 3.34 GB | `./_conda_envs/` | http://47.88.56.212/iTarget/iTarget.tar.gz |
| `python.tar.gz` | Required files for X-MOL | 1.35 GB | `./_ForFeatures/xmol/FT_to_embedding/` | http://47.88.56.212/iTarget/python.tar.gz |
| `step_400000_20200326221400.tar` | Required model files for X-MOL | 993.5 MB | `./_ForFeatures/xmol/FT_to_embedding/data/model/step_400000/` | http://47.88.56.212/iTarget/step_400000_20200326221400.tar |
| `esm2_t36_3B_UR50D.pt` | Pretrained ESM-2 weights | 5.28 GB | `./_ForFeatures/esm2/pretrained_esm2_models/` | http://47.88.56.212/iTarget/esm2_t36_3B_UR50D.pt |
| `esm2_fitted.mp` | Precomputed protein template image file | 128.3 MB | `./data/processed_data/` | http://47.88.56.212/iTarget/esm2_fitted.mp |
| `ChEMBL_data.csv` | Large-scale compound-protein interaction data collected from ChEMBL | 481.1 MB | `./data/_original_data/chembl_cpi/` | http://47.88.56.212/iTarget/ChEMBL_data.csv |

##### 3. Extract Archived Files

After downloading the compressed files, extract them into the specified directories.

```bash
# Extract python.tar.gz
tar -zxvf ./_ForFeatures/xmol/FT_to_embedding/python.tar.gz -C ./_ForFeatures/xmol/FT_to_embedding/

# Extract step_400000_20200326221400.tar
tar -xvf ./_ForFeatures/xmol/FT_to_embedding/data/model/step_400000/step_400000_20200326221400.tar -C ./_ForFeatures/xmol/FT_to_embedding/data/model/
```

##### 4. Set Up the Python Environment

iTarget mainly depends on the following package versions:

```text
python = 3.6.8
cuda = 11.1
torch = 1.8.1+cu111
biopython = 1.78
scikit-learn = 0.23.0
scipy = 1.5.4
pandas = 0.25.1
numpy = 1.19.2
lapjv = 1.3.1
umap-learn = 0.3.10
rdkit = 2019.09.3
```

Two setup options are provided. Please make sure Anaconda is installed in advance.

**Option 1**. Create the Environment from `environment.yml`

Use this option if you prefer to build the environment from scratch ( It may take a long time ).

```bash
conda env create -f environment.yml
conda activate iTarget
```

**Option 2**. Use the Packaged Conda Environment

A prebuilt Conda environment is provided for convenience.

First, locate your Anaconda installation path:

```bash
conda info | grep -i "base environment"
```

For example, the default location may be `~/anaconda3`.

Then extract the packaged environment into your own Conda environments directory:

```bash
mkdir ~/anaconda3/envs/iTarget
tar -zxvf ./_conda_envs/iTarget.tar.gz -C ~/anaconda3/envs/iTarget
```

Finally, verify that the environment is available:

```bash
conda info --envs
conda activate iTarget
```



## Training from Scratch

This section describes the full workflow for feature generation, template construction, and model training. If you **only want to use the provided pretrained model for inference**, you may skip this section and **go directly to the Inference Demo** section.

#### 1. Generate Protein and Compound Feature encodings with Language Models

This project uses ESM-2 for proteins and X-MOL for compounds.

##### 1.1 Preprocess data for benchmarks

Before running the script, modify `scripts/_data_preprocess.py` for customized usage. Supported dataset types include `bindingdb` (default), `human`, `biosnap`, and `chembl`.

```bash
bash Usage_1.1.sh
```

Generated files are saved in `./data/_original_data/`, including:

- `bindingdb_drugs.csv` and `bindingdb_prots.csv`
- `human_drugs.csv` and `human_prots.csv`
- `biosnap_drugs.csv` and `biosnap_prots.csv`
- `chembl_drugs.csv` and `chembl_prots.csv`

##### 1.2 Generate Feature Representations with ESM-2 and X-MOL

Before running this step, update `Usage_1.2.sh` for customized usage. The default settings are:

```bash
export PROTROOT=self
export DRUGROOT=self
```

The following input files are required:

- `./_ForFeatures/esm2/data/${PROTROOT}_prots.csv`
- `./_ForFeatures/xmol/FT_to_embedding/data/for_output/${DRUGROOT}_drugs.csv`

Then run:

```bash
bash Usage_1.2.sh
```

**************Notes**************: The required files for `PROTROOT=self` and `DRUGROOT=self` are already prepared in the respective target path. **If you want to generate feature representations for benchmark datasets** produced in **`Step 1.1`**, you **MUST** first copy the generated CSV files (`${PROTROOT}_prots.csv` and `${DRUGROOT}_drugs.csv`, which output and saved in `./data/_original_data` by **`step 1.1`**) into the corresponding feature-generation directories. 

**For example, for `bindingdb`:**

```bash
cp ./data/_original_data/bindingdb_drugs.csv ./_ForFeatures/xmol/FT_to_embedding/data/for_output
cp ./data/_original_data/bindingdb_prots.csv ./_ForFeatures/esm2/data
```

Then update `Usage_1.2.sh` accordingly:

```bash
export PROTROOT=bindingdb
export DRUGROOT=bindingdb
```

Then run:

```bash
bash Usage_1.2.sh
```

#### 2. Build Template Images (Optional)

This step has already been completed in the repository. You may directly use the prepared template files:

- `./data/processed_data/xmol_fitted.mp`
- `./data/processed_data/esm2_fitted.mp`

If you want to build your own custom `.mp` templates, follow the steps below:

##### 2.1 Prepare Feature Files for Template Construction

Move and rename the generated language-model feature files to the working directory `./data/original_data/scale/`.

Examples are provided for `PROTROOT=self` and `DRUGROOT=self` (`export SCALE_COURCE="self+self"`). You may modify the `{PROTROOT}, {DRUGROOT}, and {SCALE_COURCE}` setting in the `Usage_2.1.sh` and `Usage_2.2.sh` according to your own data source.

```bash
bash Usage_2.1.sh
```

##### 2.2 Compute Feature Distances and Build Template Configurations

```bash
bash Usage_2.2.sh
```

#### 3. Image-like Feature Transformation

##### 3.1 Construct Feature Images Based on the Templates with LM-feature files

This step operates on the feature files generated in **`Step 1.2`** based on the template files generated in **`Step 2.2`**. Example input files:

- `./data/original_data/example_all-data-merge-drug.csv`
- `./data/original_data/example_all-data-merge-prot.csv`

You may place your own feature files in this directory and update the `PROTROOT` and `DRUGROOT` arguments in `Usage_3.1.sh` accordingly.

**************Notes**************: Feature files generated in **`Step 1.2`** for benchmark datasets will be copied automatically to `./data/original_data/` by **`Step 1.2`**.

**For `bindingdb`, set:**

```bash
export PROTROOT=bindingdb
export DRUGROOT=bindingdb
```

Then directly run to construct Image-like features for **`bindingdb`**:

```bash
bash Usage_3.1.sh
```

#### 4. Train the Model and Run Cross-Validation

Complete this stage if you want to train iTarget from scratch and evaluate it with cross-validation.

##### 4.1 Prepare Benchmark Data for Cross-Validation

After completing **`Steps 1 to 3`** for the target benchmark, update the `BENCHMARK` argument in `Usage_4.1.sh` as needed. Supported benchmarks include `bindingdb` (default), `human`, `biosnap`, and `chembl`.

**Note**: this step is not required for `bindingdb`; you may proceed directly to **`Step 4.2`**.

```bash
bash Usage_4.1.sh
```

##### 4.2 Run Training and Cross-Validation

Update the `{ROOT}` argument in `Usage_4.2.sh` for the target benchmark. once the **`Steps 1 to 3`** for the target benchmark have been completed. The default setting is:

```bash
export ROOT=bindingdb
```

Then run:

```bash
bash Usage_4.2.sh
```

**Training results will be written to `./data/pretrained/`.**



## Inference Demo

The repository **provides a pretrained `final_model` trained on ChEMBL data**. By loading the model parameters from `./pretrained/iTarget_final_model/model.pth`, users can **predict compound-protein interactions** and **obtain the corresponding uncertainty scores** at the same time.

Given specific pair(s) of protein `sequence` and compound `SMILES`, **follow the steps below to run inference**.

#### 1. Write the sequence and SMILES into the Data Input File

Edit **`./data/predict_data/infer/infer_prot-drug.csv`** in the following format:

| source | drugid | smiles | protid | sequence | label |
| --- | --- | --- | --- | --- | --- |
| infer | inferdrug_1 | O=C......cc1 | inferprot_1 | MAS......TDY | 1 |
| infer | inferdrug_2 | NS(......)s1 | inferprot_2 | MSS......SMS | 1 |

#### 2. Generate Image-like Feature for the Input Proteins and Compounds

```bash
bash Usage_Infer1.sh
```

#### 3. Run Inference

```bash
bash Usage_Infer2.sh
```

Finally, **Inference results will be saved to `./data/predict_data/inference_results.csv`.**



## Citation and Disclaimer

The manuscript is currently under peer review.
