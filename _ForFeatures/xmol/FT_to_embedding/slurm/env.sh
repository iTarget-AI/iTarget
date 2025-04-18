#!/bin/bash 
set -xe

#export iplist=10.255.104.19,10.255.138.19,10.255.122.21

# pyhton
#export PATH=$PWD/python/bin/:$PATH
export BASE_PATH="$PWD"
# export PATH="${BASE_PATH}/python/bin/:$PATH"
# export PYTHONPATH="${BASE_PATH}/python/"
# export PATH="${BASE_PATH}/__python/bin/:$PATH"
# export PYTHONPATH="${BASE_PATH}/__python/"


#library
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
# export LD_LIBRARY_PATH=/home/work/cuda-9.0/lib64:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/home/work/cudnn/cudnn_v7/cuda/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH="${BASE_PATH}/nccl_2.3.5/lib/:$LD_LIBRARY_PATH"

#http_proxy
unset http_proxy
unset https_proxy
