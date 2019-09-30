#!/bin/sh
#SBATCH -p general
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 16000
#SBATCH -t 0-48:00
#SBATCH -J rand_treePL_bootstraps
#SBATCH -o rand_treePL_bootstraps_%A_%a.out
#SBATCH -e rand_treePL_bootstraps_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gustavo_bravo@fas.harvard.edu

module load treePL/1.0-fasrc01

./slurms_random

