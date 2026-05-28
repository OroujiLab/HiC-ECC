#!/bin/bash
#SBATCH --job-name=enhance_hg38
#SBATCH --output=logs/enhance_hg38_%j.out
#SBATCH --error=logs/enhance_hg38_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=12
#SBATCH -p veryhimem

source ~/miniconda3/etc/profile.d/conda.sh
conda activate hic

cd /cluster/home/t111631uhn/HiC-ECC
python run.py --config config/config_hg38.yaml --module enhance --sample-idx $1



#for i in {0..28}; do sbatch submit_enhance_hg38.sh $i; done

