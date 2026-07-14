#!/usr/bin/bash
set -euo pipefail

ENV_NAME="iTarget"
ENV_ARCHIVE="./_conda_envs/iTarget.tar.gz"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: Conda was not found. Please install Anaconda or Miniconda first."
    exit 1
fi

if conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
    echo "Conda environment '${ENV_NAME}' already exists."
    echo "Environment setup finished."
    exit 0
fi

CONDA_BASE="$(conda info --base)"
ENV_DIR="${CONDA_BASE}/envs/${ENV_NAME}"

if [ -f "${ENV_ARCHIVE}" ]; then
    echo "Creating Conda environment '${ENV_NAME}' from ${ENV_ARCHIVE}."
    mkdir -p "${ENV_DIR}"
    tar -zxvf "${ENV_ARCHIVE}" -C "${ENV_DIR}"
    echo "Environment setup finished."
    echo "The environment keeps the legacy name: ${ENV_NAME}."
    exit 0
fi

echo "Packaged environment archive was not found: ${ENV_ARCHIVE}"
echo "Please create the environment manually:"
echo "conda env create -f environment.yml"
echo "conda activate ${ENV_NAME}"
exit 1
