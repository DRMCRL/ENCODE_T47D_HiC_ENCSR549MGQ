import pandas as pd
import os
import re

configfile: "config/config.yml"

# Samples
df = pd.read_table(config["samples"])
samples = list(set(df['sample']))
suffix = config['suffix']

#################################
## Variables for the reference ##
#################################
ref_root = os.path.join(config['ref']['root'], "gencode-release-" + str(config['ref']['gencode']),
                        config['ref']['build'], "dna")
# Key output files
assembly = config['ref']['assembly'] + ".genome"
ref_fa = config['ref']['build'] + "." + assembly + ".fa"
ref_fagz = ref_fa + ".gz"
chr_sizes = os.path.join("output", config['ref']['build'] + ".chr_sizes.tsv")
rs_frags = os.path.join("output", config['ref']['build'] + "_" + config['hicpro']['enzyme'] + "_fragment.bed")

#####################
## HiC-Pro outputs ##
#####################
bins = re.split(r" ", config['hicpro']['bin_size'])
hicpro_config = "output/hicpro-config.txt"
HIC_PAIRS = expand(["data/hic/hic_results/data/{sample}/{sample}_allValidPairs"],
                   sample = samples)
HIC_MAT = expand(["data/hic/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}.matrix"],
                 bin = bins, sample = samples)
HIC_BED = expand(["data/hic/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}_abs.bed"],
                 bin = bins, sample = samples)
MAX_HIC = 'output/maxhic.txt'

## Define all the required outputs as a single object
REFS = [chr_sizes, rs_frags]
FAGZ = [os.path.join(ref_root, ref_fagz)]
BOWTIEIDX = expand(["{pre}.{suffix}.bt2"],
                 pre = os.path.join(ref_root, "bt2", config['ref']['build'] + "." + assembly),
                 suffix = ['1', '2', '3', '4', 'rev.1', 'rev.2'])
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
ALL_OUTPUTS.extend(BOWTIEIDX)
ALL_OUTPUTS.extend(FAGZ)
ALL_OUTPUTS.extend(FQC_OUTS)
ALL_OUTPUTS.extend(TRIM_OUTS)
ALL_OUTPUTS.extend([hicpro_config])
ALL_OUTPUTS.extend(HIC_PAIRS)
ALL_OUTPUTS.extend(HIC_MAT)
ALL_OUTPUTS.extend(HIC_BED)
ALL_OUTPUTS.extend([MAX_HIC])

rule all:
    input:
        ALL_OUTPUTS

include: "rules/download.smk"
include: "rules/reference.smk"
include: "rules/qc.smk"
include: "rules/trimming.smk"
include: "rules/hicpro.smk"
include: "rules/maxhic.smk"

