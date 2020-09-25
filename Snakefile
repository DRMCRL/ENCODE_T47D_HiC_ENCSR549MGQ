import pandas as pd

configfile: "config/config.yml"

# Samples
samples = pd.read_table(config["samples"])

## Variables for the reference
ref_root = config['ref']['root'] + "/gencode-release-" + config['ref']['gencode'] + "/" + config['ref']['build'] + "/dna"
ref_fagz = config['ref']['build'] + ".primary_assembly.genome.fa.gz"
chr_sizes = ref_root + "/" + config['ref']['build'] + ".chr_sizes.tsv"
rs_frags = ref_root + "/" + config['ref']['build'] + "_" + config['ref']['enzyme'] + "_fragment.bed"

# Define all the required outputs as a single object
REF = ref_root + "/" + ref_fagz
BOWTIEIDX = expand(["{path}/{build}.primary_assembly.genome.{suffix}.bt2"],
                   path = ref_index,
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
ALL_OUTPUTS = []
ALL_OUTPUTS.extend(REF)
ALL_OUTPUTS.extend(BOWTIEIDX)
ALL_OUTPUTS.extend(chr_sizes)
ALL_OUTPUTS.extend(rs_frags)
ALL_OUTPUTS.extend(FQC_OUTS)
ALL_OUTPUTS.extend(TRIM_OUTS)

rule all:
    input:
        ALL_OUTPUTS

include: "reference.smk"
include: "rules/qc.smk"
include: "rules/trimming.smk"


