rule get_chrom_sizes:
    output: chr_sizes
    params:
        ftp = "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/all_assembly_versions",
        genbank = config['ref']['genbank'],
        build = config['ref']['build']
    threads: 1
    shell:
        """

        # Download the assembly report
        TEMPDIR=$(mktemp -d -t chrXXXXXXXXXX)
        REPORT="{params.genbank}_{params.build}_assembly_report.txt"
        URL="{params.ftp}/{params.genbank}_{params.build}/{params.genbank}_{params.build}_assembly_report.txt"
        wget -O "$TEMPDIR/$REPORT" $URL

        # Extract the chrom_sizes
        egrep 'assembled-molecule' "$TEMPDIR/$REPORT" | \
          awk '{{print $11"\t"$10}}' > {output}

        rm -rf $TEMPDIR

        """

rule find_rs_fragments:
    input: os.path.join(ref_root, ref_fa)
    output: rs_frags
    params:
        enzyme = config['hicpro']['enzyme']
    threads: 1
    conda: "../envs/python2.7.yml"
    shell:
        """
        # Get the latest version from the HiC-Pro repo
        wget \
          -O "scripts/digest_genome.py" \
          "https://raw.githubusercontent.com/nservant/HiC-Pro/master/bin/utils/digest_genome.py"

        # Run the python script
        python scripts/digest_genome.py \
          -r {params.enzyme} \
          -o {output} \
          {input}
        """

rule make_hicpro_config:
    input:
        idx = os.path.join(ref_root, "bt2"),
        rs = rs_frags,
        chr_sizes = chr_sizes
    output:
        hicpro_config
    conda: "../envs/tidyverse.yml"
    threads: 1
    shell:
        """
        Rscript --vanilla \
          scripts/write_hicpro_config.R \
          {input.idx} \
          {input.chr_sizes} \
          {input.rs} \
          {output}
        """

rule run_hicpro:
    input:
        config = hicpro_config,
        files = expand(["data/trimmed/fastq/{sample}/{sample}{reads}{suffix}"],
                       sample = samples, suffix = suffix,
                       reads = [config['hicpro']['pair1_ext'], config['hicpro']['pair2_ext']])
    output:
        valid_pairs = expand(["data/hic/hic_results/data/{sample}/{sample}_allValidPairs"],
                             sample = samples),
        mat = expand(["data/hic/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}.matrix"],
                     sample = samples, bin = bins),
        bed = expand(["data/hic/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}_abs.bed"],
                     sample = samples, bin = bins)
    log: "logs/hicpro/run_hicpro.log"
    threads: config['hicpro']['ncpu']
    shell:
        """
        ######################################
        ## Specific to phoenix for now only ##
        ######################################
        ## Load modules
        module load HiC-Pro/2.9.0-foss-2016b

        ##Run HiC-pro
        HiC-Pro \
          -c {input.config} \
          -i "data/trimmed/fastq" \
          -o "data/hic" &> {log}
        """



