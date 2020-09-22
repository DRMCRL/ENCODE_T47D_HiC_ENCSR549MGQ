#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=4:00:00
#SBATCH --mem=64GB
#SBATCH -o /home/a1018048/slurm/ENCODE_T47D_HiC_ENCSR549MGQ/%x_%j.out
#SBATCH -e /home/a1018048/slurm/ENCODE_T47D_HiC_ENCSR549MGQ/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=a1018048@adelaide.edu.au

# to prepare all the reference files for HiC-Pro
module load SAMtools/1.3.1-foss-2016b
module load Bowtie2/2.2.9-foss-2016b
module load Python/2.7.13-foss-2016b

# Define the directories
USR=/hpcfs/users/a1018048
PROJROOT=${USR}/ENCODE_T47D_HiC_ENCSR549MGQ

# Setup the assembly names & details
BUILD=GRCh37
ID=GCA_000001405.1
ASSEMB=${ID}_${BUILD}
SP=Homo_sapiens
REPORT=${ASSEMB}_assembly_report.txt
FTP=ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/${SP}/all_assembly_versions

# Check then Unzip the reference
REFBASE=${USR}/refs/genocode-release-33/${BUILD}/dna
GZREF=${REFBASE}/${BUILD}.primary_assembly.genome.fa.gz
FAREF=${GZREF%.gz}
if [[ ! -f ${GZREF} ]]; then
  echo "Couldn't find ${GZREF}"
  exit 1
fi
if [[ -f ${FAREF} ]]; then
  echo "Found ${FAREF}\n"
else
  echo "Unzipping ${GZREF}"
  gunzip -c ${GZREF} > ${FAREF}
fi

# build the bowtie2 index
mkdir -p ${REFBASE}/bt2
bowtie2-build \
  -f ${FAREF} \
  --threads 16 \
  ${REFBASE}/bt2/$(basename ${FAREF%.fa})

# get chrom size
wget \
  -O ${REFBASE}/${REPORT} \
  ${FTP}/${ASSEMB}/${REPORT}

egrep 'assembled-molecule' ${REPORT} | awk '{print $11"\t"$10}' > ${BUILD}.chr_sizes.tsv

# get restriction fragments
python ${PROJROOT}/scripts/digest_genome.py \
  -r A^AGCTT \
  -o ${REFBASE}/GRCh37_HindIII_fragment.bed \
  ${FAREF}

# Remove the unzippped reference
rm ${FAREF}
