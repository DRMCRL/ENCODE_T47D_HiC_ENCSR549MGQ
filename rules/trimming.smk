rule adapter_removal:
    input:
        r1 = "data/raw/fastq/{sample}_R1.fastq.gz",
        r2 = "data/raw/fastq/{sample}_R2.fastq.gz"
    output:
        t1 = "data/trimmed/fastq/{sample}_R1.fastq.gz",
        t2 = "data/trimmed/fastq/{sample}_R2.fastq.gz",
        log = "data/trimmed/logs/{sample}.settings"
    conda:
        "../envs/adapterremoval.yml"
    params:
        adapter1 = config['trimming']['adapter1'],
        adapter2 = config['trimming']['adapter2'],
        minlength = config['trimming']['minlength'],
        minqual = config['trimming']['minqual'],
        maxns = config['trimming']['maxns']
    threads: 1
    log:
        "logs/adapterremoval/{sample}.log"
    shell:
        """
        SAMPLE=$(basename {input.r1})
        DISCARDDIR="data/trimmed/discarded"
        DISCARD="${DISCARDDIR}/${SAMPLE%_R1.fastq.gz}.discarded.fastq.gz"
        SINGLE="${DISCARDDIR}/${SAMPLE%_R1.fastq.gz}.singleton.truncated.fastq.gz"
        AdapterRemoval \
            --adapter1 {params.adapter1} \
            --adapter2 {params.adapter2} \
            --file1 {input.r1} \
            --file2 {input.r2} \
            --threads {threads} \
            --gzip \
            --maxns {params.maxns} \
            --trimqualities \
            --minquality {params.minqual} \
            --minlength {params.minlength} \
            --output1 {output.t1} \
            --output2 {output.t2} \
            --discarded ${DISCARD} \
            --singleton ${SINGLE} \
            --settings {output.log} &> {log}
        """
