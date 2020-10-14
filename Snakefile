import pandas as pd
import os
import re

configfile: "config/config.yml"

# Samples
samples = pd.read_table(config["samples"])

## Variables for the reference
ref_root = os.path.join(config['ref']['root'], "gencode-release-" + str(config['ref']['gencode']),
                        config['ref']['build'], "dna")
# Key output files
ref_fa = config['ref']['build'] + ".primary_assembly.genome.fa"
ref_fagz = ref_fa + ".gz"
chr_sizes = os.path.join(ref_root, config['ref']['build'] + ".chr_sizes.tsv")
rs_frags = os.path.join(ref_root, config['ref']['build'] + "_" + config['hicpro']['enzyme'] + "_fragment.bed")

## Define all the required outputs as a single object
REFS = expand(["{ref_root}/{build}{suffix}"],
              ref_root = ref_root,
              build = config['ref']['build'],
              suffix = ['.chr_sizes.tsv', "_" + config['hicpro']['enzyme'] + "_fragment.bed"])
FAGZ = expand(["{path}/{file}"], path = ref_root, file = ref_fagz)
BOWTIEIDX = expand(["{path}/{build}.primary_assembly.genome.{suffix}.bt2"],
                   path = os.path.join(ref_root, "bt2"),
                   build = config['ref']['build'],
                   suffix = ['1', '2', '3', '4', 'rev.1', 'rev.2'])
FQC_OUTS = expand(["data/{step}/FastQC/{sample}_{reads}_fastqc.{suffix}"],
                 suffix = ['zip', 'html'],
                 reads = ['R1', 'R2'],
                 sample = samples['sample'],
                 step = ['raw', 'trimmed'])
TRIM_OUTS = expand(["data/trimmed/fastq/{sample}/{sample}_{reads}.fastq.gz"],
                  sample = samples['sample'],
                  reads = ['R1', 'R2'])
HIC_CONFIG = ['config/hicpro-config.txt']
ALL_OUTPUTS = []
ALL_OUTPUTS.extend(REFS)
ALL_OUTPUTS.extend(BOWTIEIDX)
ALL_OUTPUTS.extend(FAGZ)
ALL_OUTPUTS.extend(FQC_OUTS)
ALL_OUTPUTS.extend(TRIM_OUTS)
ALL_OUTPUTS.extend(HIC_CONFIG)

rule all:
    input:
        ALL_OUTPUTS

include: "rules/reference.smk"
include: "rules/qc.smk"
include: "rules/trimming.smk"
include: "rules/hicpro.smk"


