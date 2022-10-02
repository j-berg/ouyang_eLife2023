#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH -o /uufs/chpc.utah.edu/common/home/rutter-group2/j-berg/ahmad_19175R/slurmjob-%j
#SBATCH --account=rutter-gpu-np
#SBATCH --partition=rutter-gpu-np
#SBATCH --gres=gpu:2080ti:4
#SBATCH --mem=0

# Rutter lab GPU node specs 
# - 40 cores
# - 192 GB of memory
# - 4x RTX2080TI GPUs

# Note: Change --gres command to up to `:4` to gain access to all 4 GPUs on node
# Note: --mem=0 reserves all node memory for the current job 
# See https://www.chpc.utah.edu/documentation/guides/gpus-accelerators.php for more information on options
echo "+ Checking CUDA devices..."
echo "Devices: $CUDA_VISIBLE_DEVICES"
nvidia-smi -L

NCPU=$(grep -c ^processor /proc/cpuinfo)
echo "No. CPUs: $NCPU"

cd /uufs/chpc.utah.edu/common/home/rutter-group2/j-berg/yeyun_2022_04_13_novogene
wget -r -c ftp://X202SC22032816-Z01-F001:k9jt99a8@usftp21.novogene.com:21/
