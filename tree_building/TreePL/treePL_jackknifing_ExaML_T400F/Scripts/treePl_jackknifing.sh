#!/bin/sh
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 8000
#SBATCH -t 0-48:00
#SBATCH -J treePL_jackknifing
#SBATCH -o treePL_jackknifing_%A_%a.out
#SBATCH -e treePL_jackknifing_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gustavo_bravo@fas.harvard.edu

module load treePL/1.0-fasrc01

./slurms_jackknifing
