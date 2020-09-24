#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=48:00:00
#SBATCH --mem=64GB
#SBATCH -o /home/a1018048/slurm/ENCODE_T47D_HiC_ENCSR549MGQ/%x_%j.out
#SBATCH -e /home/a1018048/slurm/ENCODE_T47D_HiC_ENCSR549MGQ/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stephen.pederson@adelaide.edu.au

## Bash script to run HiC-pro
## By Ning Liu, modified by Steve Pederson

## Load modules
module load HiC-Pro/2.9.0-foss-2016b

## Define the root for the analysis
ROOTDIR=/hpcfs/users/a1018048/ENCODE_T47D_HiC_ENCSR549MGQ
DATADIR=${ROOTDIR}/data/trimmed/fastq
HICDIR=${ROOTDIR}/data/hic
CONFIG=${ROOTDIR}/config/hicpro-config.txt

##Run HiC-pro
HiC-Pro \
  -c ${CONFIG} \
  -i ${DATADIR} \
  -o ${HICDIR}
