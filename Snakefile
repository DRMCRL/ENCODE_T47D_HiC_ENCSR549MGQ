import pandas as pd

configfile: "config/config.yml"

# Samples
samples = pd.read_table(config["samples"])

## Variables for the reference
ref_root = config['ref']['root']
ref_index = ref_root + "/" + config['ref']['index']

# Define all the required outputs as a single object
FQC_OUTS = expand(["data/{step}/FastQC/{sample}_fastqc.{suffix}"],
                 sample = samples['file'],
                 suffix = ["zip", "html"],
                 step = ["raw", "trimmed"])
TRIM_OUTS = expand(["data/trimmed/fastq/{sample}_{reads}.fastq.gz"],
                  sample = samples['sample'],
                  reads = ["R1", "R2"])
ALL_OUTPUTS = []
ALL_OUTPUTS.extend(FQC_OUTS)
ALL_OUTPUTS.extend(TRIM_OUTS)

rule all:
    input:
        ALL_OUTPUTS

include: "rules/qc.smk"
include: "rules/trimming.smk"


