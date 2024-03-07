#! /bin/bash

#SBATCH --job-name=test
#SBATCH --partition=nice
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=test_%j.log

module load python/3.10.12-qdw4iom
module load foldseek/8-ef4e960-5kdayhs
module load miniconda3
conda activate pyFoldSeek
