import pandas as pd
import os
import re

configfile: "config/config.yml"

# Samples
df = pd.read_table(config["samples"])
samples = list(set(df['sample']))
suffix = config['suffix']
read_ext = [config['hicpro']['pair1_ext'], config['hicpro']['pair2_ext']]

#################################
## Variables for the reference ##
#################################

build = config['ref']['build']
ref_root = os.path.join(config['ref']['root'], "gencode-release-" + str(config['ref']['gencode']),
                        build, "dna")
# Key output files
assembly = config['ref']['assembly'] + ".genome"
ref_fa = build + "." + assembly + ".fa"
ref_fagz = ref_fa + ".gz"
chr_sizes = os.path.join(os.getcwd(), "config", build + ".chr_sizes.tsv")
rs_frags = os.path.join(os.getcwd(), "config", build + "_" + config['hicpro']['enzyme'] + "_fragment.bed")

##########################################################
## Define all the required outputs from the setup steps ##
##########################################################

REFS = [chr_sizes, rs_frags]
FAGZ = [os.path.join(ref_root, ref_fagz)]
BOWTIEIDX = expand([ref_root + "/bt2/{prefix}.{sub}.bt2"],
               prefix = config['ref']['build'] + "." + assembly,
               sub = ['1', '2', '3', '4', 'rev.1', 'rev.2'] )
FQC_OUTS = expand(["data/{step}/FastQC/{sample}_{reads}_fastqc.{suffix}"],
                 suffix = ['zip', 'html'],
                 reads = ['R1', 'R2'],
                 sample = samples,
                 step = ['raw', 'trimmed'])
TRIM_OUTS = expand(["data/trimmed/fastq/{sample}/{sample}_{reads}{suffix}"],
                  sample = samples, suffix = suffix,
                  reads = ['R1', 'R2'])

ALL_OUTPUTS = []
ALL_OUTPUTS.extend(REFS)
ALL_OUTPUTS.extend([BOWTIEIDX])
ALL_OUTPUTS.extend(FAGZ)
ALL_OUTPUTS.extend(FQC_OUTS)
ALL_OUTPUTS.extend(TRIM_OUTS)

#####################
## HiC-Pro outputs ##
#####################
bins = re.split(r" ", config['hicpro']['bin_size'])
hicpro_config = "config/hicpro-config.txt"
digest_script = "scripts/digest_genome.py"
MAPPING =  expand(["data/hic/bowtie_results/bwt2/{sample}/{sample}{reads}_" + build + "." + assembly + ".bwt2merged.bam"],
                  reads = read_ext, sample = samples)
PROC_BAM = expand(["data/hic/bowtie_results/bwt2/{sample}/{sample}_" + build + "." + assembly + ".bwt2pairs.bam"],
                  sample = samples)
PROC_PAIRS = expand(["data/hic/hic_results/data/{sample}/{sample}_" + build + "." + assembly + ".bwt2pairs.validPairs"],
                    sample = samples)
HIC_QC = ['data/hic/hic_results/pic']
VALID_PAIRS = expand(["data/hic/hic_results/data/{sample}/{sample}_allValidPairs"],
                       sample = samples)
HIC_MAT = expand(["data/hic/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}.matrix"],
                  sample = samples, bin = bins)
HIC_BED = expand(["data/hic/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}_{type}.bed"],
                  sample = samples, bin = bins, type = ['abs', 'ord'])

ALL_OUTPUTS.extend([hicpro_config, digest_script])
ALL_OUTPUTS.extend(MAPPING)
ALL_OUTPUTS.extend(PROC_BAM)
ALL_OUTPUTS.extend(PROC_PAIRS)
ALL_OUTPUTS.extend(HIC_QC)
ALL_OUTPUTS.extend(VALID_PAIRS)
ALL_OUTPUTS.extend(HIC_MAT)
ALL_OUTPUTS.extend(HIC_BED)

#####################
## Max HiC Outputs ##
#####################
MAXHIC_INTERACTIONS = expand(["output/MaxHiC/{sample}/{bin}/{type}_interactions.txt"],
                             sample = samples, bin = ['40000'],
                             type = ['cis', 'trans'])
ALL_OUTPUTS.extend(MAXHIC_INTERACTIONS)

#####################
## Rules & Outputs ##
#####################

rule all:
    input:
        ALL_OUTPUTS

include: "rules/download.smk"
include: "rules/reference.smk"
include: "rules/qc.smk"
include: "rules/trimming.smk"
include: "rules/hicpro.smk"
include: "rules/maxhic.smk"

