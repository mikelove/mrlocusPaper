#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=240
#SBATCH --mem=1000

module load anaconda/2019.10
module load r/4.0.1
module load gcc/6.3.0

snakemake -j 10 --latency-wait 30 --cluster "sbatch --mem=5000 -t 5"
