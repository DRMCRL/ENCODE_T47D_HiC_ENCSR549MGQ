rule raw_fastqc:
    input:
        "data/raw/fastq/{sample}.fastq.gz"
    output:
        html = "data/raw/FastQC/{sample}_fastqc.html",
        zip = "data/raw/FastQC/{sample}_fastqc.zip"
    conda:
        "../envs/fastqc.yml"
    params: config['fastqc']['params']
    log:
        "logs/FastQC/raw/{sample}.log"
    threads: 1
    wrapper:
        "0.66.0/bio/fastqc"

rule trim_fastqc:
    input:
        "data/trimmed/fastq/{sample}.fastq.gz"
    output:
        html = "data/trimmed/FastQC/{sample}_fastqc.html",
        zip = "data/trimmed/FastQC/{sample}_fastqc.zip"
    conda:
        "../envs/fastqc.yml"
    params: config['fastqc']['params']
    log:
        "logs/FastQC/trimmed/{sample}.log"
    threads: 1
    wrapper:
        "0.66.0/bio/fastqc"
