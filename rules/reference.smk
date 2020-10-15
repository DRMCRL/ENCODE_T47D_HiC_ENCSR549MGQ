rule get_reference:
    output: os.path.join(ref_root, ref_fa)
    params:
        genbank = config['ref']['genbank'],
        gencode = config['ref']['gencode'],
        build = config['ref']['build']
    threads: 1
    shell:
        """
        # Define the URL and download
        URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{params.gencode}/{params.build}_mapping/$(basename {output}).gz"
        wget \
          -O "{output}.gz" \
          $URL
        gunzip -c "{output}.gz" > {output}
        """

rule bowtie2_index:
    input: os.path.join(ref_root, ref_fa)
    output:
        dir = os.path.join(ref_root, "bt2")
    conda: "../envs/bowtie2.yml"
    threads: 8
    params:
        prefix = config['ref']['build'] + "." + assembly
    shell:
        """
        bowtie2-build \
          --threads {threads} \
          -f {input} \
          {output}/{params.prefix}
        """

rule rezip_fa:
    input:
        frags = rs_frags,
        temp_fa = rules.get_reference.output,
        dir = os.path.join(ref_root, "bt2")
    output: os.path.join(ref_root, ref_fagz)
    shell:
        """
        gzip -c {input.temp_fa} > {output}
        rm {input.temp_fa}
        """
