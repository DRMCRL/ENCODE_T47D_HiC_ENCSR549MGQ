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
    params:
        ncpu = config['hicpro']['ncpu'],
        pair1_ext = config['hicpro']['pair1_ext'],
        pair2_ext = config['hicpro']['pair2_ext'],
        min_mapq = config['hicpro']['min_mapq'],
        phred = config['hicpro']['phred'],
        genome = config['ref']['build'] + "." + assembly,
        ligation_site = config['hicpro']['ligation_site'],
        bin_size = config['hicpro']['bin_size']
    threads: 1
    shell:
        """
        echo -e "## Paths and Settings" >> {output}
        echo -e "TMP_DIR = tmp\nLOGS_DIR = logs\nBOWTIE2_OUTPUT_DIR = bowtie_results\nMAPC_OUTPUT = hic_results\nRAW_DIR = rawdata" >> {output}

        echo -e "\n##SYSTEM AND SCHEDULER" >> {output}
        echo -e "N_CPU = {params.ncpu}" >> {output}
        echo -e "LOGILE = hicpro.log" >> {output}
        echo -e "JOB_NAME =\nJOB_MEM =\nJOB_WALLTIME =\nJOB_QUEUE =\nJOB_MAIL =" >> {output}

        echo -e "\n## Data" >> {output}
        echo -e "PAIR1_EXT = {params.pair1_ext}" >> {output}
        echo -e "PAIR2_EXT = {params.pair2_ext}" >> {output}

        echo -e "\n## Alignment Options" >> {output}
        echo -e "FORMAT = phred{params.phred}" >> {output}
        echo -e "MIN_MAPQ = {params.min_mapq}" >> {output}
        echo -e "BOWTIE2_IDX_PATH = {input.idx}" >> {output}
        echo -e "BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder" >> {output}
        echo -e "BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder" >> {output}

        echo -e "\n## Annotation files" >> {output}
        echo -e "REFERENCE_GENOME = {params.genome}" >> {output}
        echo -e "GENOME_SIZE = {input.chr_sizes}" >> {output}
        echo -e "CAPTURE_TARGET =" >> {output}

        echo -e "\n## Allele specific analysis" >> {output}
        echo -e "ALLELE_SPECIFIC_SNP =" >> {output}

        echo -e "\n## Digestion Hi-C" >> {output}
        echo -e "GENOME_FRAGMENT = {input.rs}" >> {output}
        echo -e "LIGATION_SITE = {params.ligation_site}" >> {output}
        echo -e "MIN_FRAG_SIZE =\nMAX_FRAG_SIZE =\nMIN_INSERT_SIZE =\nMAX_INSERT_SIZE =" >> {output}

        echo -e "\n## Hi-C processing" >> {output}
        echo -e "MIN_CIS_DIST =" >> {output}
        echo -e "GET_ALL_INTERACTION_CLASSES = 1" >> {output}
        echo -e "GET_PROCESS_SAM = 0" >> {output}
        echo -e "RM_SINGLETON = 1" >> {output}
        echo -e "RM_MULTI = 1" >> {output}
        echo -e "RM_DUP = 1" >> {output}

        echo -e "\n## Contact Maps" >> {output}
        echo -e "BIN_SIZE = {params.bin_size}" >> {output}
        echo -e "MATRIX_FORMAT = complete" >> {output}

        echo -e "\n## Normalization" >> {output}
        echo -e "MAX_ITER = 100" >> {output}
        echo -e "FILTER_LOW_COUNT_PERC = 0.02" >> {output}
        echo -e "FILTER_HIGH_COUNT_PERC = 0" >> {output}
        echo -e "EPS = 0.1" >> {output}

        """

rule run_hicpro:
    input:
        config = hicpro_config,
        files = expand(["data/trimmed/fastq/{sample}/{sample}{reads}{suffix}"],
                       sample = samples, suffix = suffix,
                       reads = [config['hicpro']['pair1_ext'], config['hicpro']['pair1_ext']])
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

        ##Run HiC-pro responding to yes to any interactive requests
        yes | HiC-Pro \
          -c {input.config} \
          -i "data/trimmed/fastq" \
          -o "data/hic" &> {log}
        """



