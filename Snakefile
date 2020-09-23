import pandas as pd

configfile: "config/config.yml"

# Samples
samples = pd.read_table(config["samples"])

## Variables for the reference
ref_root = config['ref']['root']
ref_index = ref_root + "/" + config['ref']['index']

# Define all the required outputs as a single object
FQC_OUTS = expand(["data/{step}/FastQC/{sample}/{file}_fastqc.{suffix}"],
                 suffix = ["zip", "html"],
                 file = samples['file'],
                 sample = samples['sample'],
                 step = ["raw", "trimmed"])
TRIM_OUTS = expand(["data/trimmed/fastq/{sample}/{file}.fastq.gz"],
                  sample = samples['sample'],
                  file = samples['file'])
ALL_OUTPUTS = []
ALL_OUTPUTS.extend(FQC_OUTS)
ALL_OUTPUTS.extend(TRIM_OUTS)

rule all:
    input:
        ALL_OUTPUTS

include: "rules/qc.smk"
include: "rules/trimming.smk"


