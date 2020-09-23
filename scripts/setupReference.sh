#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2:00:00
#SBATCH --mem=64GB
#SBATCH -o /home/a1018048/slurm/ENCODE_T47D_HiC_ENCSR549MGQ/%x_%j.out
#SBATCH -e /home/a1018048/slurm/ENCODE_T47D_HiC_ENCSR549MGQ/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stephen.pederson@adelaide.edu.au

# to prepare all the reference files for HiC-Pro
module load SAMtools/1.3.1-foss-2016b
module load Bowtie2/2.2.9-foss-2016b
module load Python/2.7.13-foss-2016b

CORES=8

# Define the directories
USR=/hpcfs/users/a1018048
PROJROOT=${USR}/ENCODE_T47D_HiC_ENCSR549MGQ

# Define the assembly names & details
BUILD=GRCh37
ID=GCA_000001405.1
ASSEMB=${ID}_${BUILD}
SP=Homo_sapiens
REPORT=${ASSEMB}_assembly_report.txt
NCBIFTP=ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/${SP}/all_assembly_versions

#########################
## Setup the reference ##
#########################
REL=33
GENCODEFTP=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${REL}/${BUILD}_mapping/
REFBASE=${USR}/refs/gencode-release-${REL}/${BUILD}/dna
GZREF=${REFBASE}/${BUILD}.primary_assembly.genome.fa.gz
if [[ ! -d ${REFBASE} ]]; then
  mkdir -p ${REFBASE}
  wget \
    -O ${GZREF} \
    ${GENCODEFTP}/$(basename ${GZREF})
fi
if [[ ! -f ${GZREF} ]]; then
  echo "Couldn't find ${GZREF}"
  printenv
  exit 1
fi

# Unzip the reference
FAREF=${GZREF%.gz}
if [[ -f ${FAREF} ]]; then
  echo "Found ${FAREF}\n"
else
  echo "Unzipping ${GZREF}"
  gunzip -c ${GZREF} > ${FAREF}
fi

#############################
## Build the bowtie2 index ##
#############################
if [[ ! -d ${REFBASE}/bt2 ]]; then
  mkdir -p ${REFBASE}/bt2
fi
bowtie2-build \
  -f ${FAREF} \
  --threads ${CORES} \
  ${REFBASE}/bt2/$(basename ${FAREF%.fa})

# get chrom size
wget \
  -O ${REFBASE}/${REPORT} \
  ${NCBIFTP}/${ASSEMB}/${REPORT}

egrep 'assembled-molecule' ${REFBASE}/${REPORT} | \
  awk '{print $11"\t"$10}' > \
  ${REFBASE}/${BUILD}.chr_sizes.tsv

###############################
## Get restriction fragments ##
###############################
RSFILE=${REFBASE}/${BUILD}_HindIII_fragment.bed
python ${PROJROOT}/scripts/digest_genome.py \
  -r A^AGCTT \
  -o ${RSFILE} \
  ${FAREF}

# Remove the unzippped reference
rm ${FAREF}
