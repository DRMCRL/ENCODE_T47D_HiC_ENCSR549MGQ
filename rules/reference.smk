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
    input: rules.get_reference.output
    output:
        bt2 = expand(["{pre}.{suffix}.bt2"],
                 pre = os.path.join(ref_root, "bt2", config['ref']['build'] + "." + assembly),
                 suffix = ['1', '2', '3', '4', 'rev.1', 'rev.2'])
    conda: "../envs/bowtie2.yml"
    threads: 8
    params:
        idx_root = os.path.join(ref_root, "bt2"),
        prefix = config['ref']['build'] + "." + assembly
    shell:
        """
        bowtie2-build \
          --threads {threads} \
          -f {input} \
          {params.idx_root}/{params.prefix}
        """

rule rezip_fa:
    input:
        frags = rs_frags,
        temp_fa = rules.get_reference.output,
        bt2 = rules.bowtie2_index.output.bt2
    output: os.path.join(ref_root, ref_fagz)
    shell:
        """
        gzip -c {input.temp_fa} > {output}
        rm {input.temp_fa}
        """
