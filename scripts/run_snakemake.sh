#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=1:30:00
#SBATCH --mem=32GB
#SBATCH -o /home/a1018048/slurm/ENCODE_T47D_HiC_ENCSR549MGQ/%x_%j.out
#SBATCH -e /home/a1018048/slurm/ENCODE_T47D_HiC_ENCSR549MGQ/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stephen.pederson@adelaide.edu.au

## Cores
CORES=8
if [ -d "/hpcfs" ]; then
	module load arch/arch/haswell
	module load arch/haswell
	module load modulefiles/arch/haswell
	HPC="/hpcfs"
else
    if [ -d "/fast" ]; then
        HPC=/fast
    else
        exit 1
    fi
fi

## Project Root
PROJ=${HPC}/users/a1018048/ENCODE_T47D_HiC_ENCSR549MGQ

## The environment containing snakemake
micromamba activate snakemake
cd ${PROJ}

## Create dot and pdf files for visualisation
#snakemake --rulegraph > output/rulegraph.dot
#dot -Tpdf output/rulegraph.dot > output/rulegraph.pdf

## Run snakemake
snakemake \
  --cores ${CORES} \
  --use-conda \
  --latency-wait 300 \
  -- delete-temp-output \
  --until make_hicpro_config
