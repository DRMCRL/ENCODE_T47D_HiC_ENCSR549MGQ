rule make_hicpro_config:
    output:
        "config/hicpro-config.txt"
    params:
        ncpu = config['hicpro']['ncpu'],
        pair1_ext = config['hicpro']['pair1_ext'],
        pair2_ext = config['hicpro']['pair2_ext'],
        min_mapq = config['hicpro']['min_mapq'],
        phred = config['hicpro']['phred'],
        idx = os.path.join(ref_root, "bt2"),
        genome = config['ref']['build'] + "." + assembly,
        chr_sizes = chr_sizes,
        fragment_bed = rs_frags,
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
        echo -e "BOWTIE2_IDX_PATH = {params.idx}" >> {output}
        echo -e "BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder" >> {output}
        echo -e "BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder" >> {output}

        echo -e "\n## Annotation files" >> {output}
        echo -e "REFERENCE_GENOME = {params.genome}" >> {output}
        echo -e "GENOME_SIZE = {params.chr_sizes}" >> {output}
        echo -e "CAPTURE_TARGET =" >> {output}

        echo -e "\n## Allele specific analysis" >> {output}
        echo -e "ALLELE_SPECIFIC_SNP =" >> {output}

        echo -e "\n## Digestion Hi-C" >> {output}
        echo -e "GENOME_FRAGMENT = {params.fragment_bed}" >> {output}
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
