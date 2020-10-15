import pandas as pd
import os
import re

configfile: "config/config.yml"

# Samples
df = pd.read_table(config["samples"])
samples = list(set(df['sample']))

## Variables for the reference
ref_root = os.path.join(config['ref']['root'], "gencode-release-" + str(config['ref']['gencode']),
                        config['ref']['build'], "dna")
# Key output files
assembly = config['ref']['assembly'] + ".genome"
ref_fa = config['ref']['build'] + "." + assembly + ".fa"
ref_fagz = ref_fa + ".gz"
chr_sizes = os.path.join(ref_root, config['ref']['build'] + ".chr_sizes.tsv")
rs_frags = os.path.join(ref_root, config['ref']['build'] + "_" + config['hicpro']['enzyme'] + "_fragment.bed")

## HiC-Pro outputs
bins = re.split(r" ", config['hicpro']['bin_size'])
hic_dir = "data/hic"
HIC_CONFIG = ['config/hicpro-config.txt']
HIC_PAIRS = expand(["{path}/hic_results/data/{sample}/{sample}_allValidPairs"],
                   sample = samples, path = hic_dir)
HIC_MAT = expand([os.path.join(hic_dir, "hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}{suffix}")],
                      suffix = ['_abs.bed', '.matrix'],
                      bin = bins,
                      sample = samples)

## Define all the required outputs as a single object
REFS = expand(["{ref_root}/{build}{suffix}"],
              ref_root = ref_root,
              build = config['ref']['build'],
              suffix = ['.chr_sizes.tsv', "_" + config['hicpro']['enzyme'] + "_fragment.bed"])
FAGZ = [os.path.join(ref_root, ref_fagz)]
BOWTIEIDX = expand(["{pre}.{suffix}.bt2"],
                 pre = os.path.join(ref_root, "bt2", config['ref']['build'] + "." + assembly),
                 suffix = ['1', '2', '3', '4', 'rev.1', 'rev.2'])
FQC_OUTS = expand(["data/{step}/FastQC/{sample}_{reads}_fastqc.{suffix}"],
                 suffix = ['zip', 'html'],
                 reads = ['R1', 'R2'],
                 sample = samples,
                 step = ['raw', 'trimmed'])
TRIM_OUTS = expand(["data/trimmed/fastq/{sample}/{sample}_{reads}.fastq.gz"],
                  sample = samples,
                  reads = ['R1', 'R2'])
ALL_OUTPUTS = []
ALL_OUTPUTS.extend(REFS)
ALL_OUTPUTS.extend(BOWTIEIDX)
ALL_OUTPUTS.extend(FAGZ)
ALL_OUTPUTS.extend(FQC_OUTS)
ALL_OUTPUTS.extend(TRIM_OUTS)
ALL_OUTPUTS.extend(HIC_CONFIG)
ALL_OUTPUTS.extend(HIC_PAIRS)
ALL_OUTPUTS.extend(HIC_MAT)
ALL_OUTPUTS.extend([hic_dir])

rule all:
    input:
        ALL_OUTPUTS

include: "rules/reference.smk"
include: "rules/qc.smk"
include: "rules/trimming.smk"
include: "rules/hicpro.smk"


