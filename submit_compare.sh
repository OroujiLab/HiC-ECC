#!/bin/bash
#SBATCH --job-name=chess_compare
#SBATCH --output=logs/chess_%j.out
#SBATCH --error=logs/chess_%j.err
#SBATCH --time=36:00:00
#SBATCH --mem=90G
#SBATCH --cpus-per-task=4
#SBATCH -p veryhimem

#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate hic

cd /cluster/home/t111631uhn/HiC-ECC
#python run.py --config config/config.yaml --module compare --tissue-idx $1

python run.py --config config/config.yaml --module enhance



#for i in {0..7}; do sbatch submit_compare.sh $i; done

