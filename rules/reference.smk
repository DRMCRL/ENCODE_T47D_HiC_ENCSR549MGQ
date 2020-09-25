rule get_reference:
    output: ref_root + "/" + ref_fagz
    params:
        genbank = config['ref']['genbank'],
        gencode = config['ref']['gencode'],
        build = config['ref']['build']
    threads: 1
    shell:
        """
        # Define the URL and download
        URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{params.gencode}/{params.build}_mapping/$(basename {output})"
        wget \
          -O {output} \
          $URL
        """

rule bowtie2_index:
    input: rules.get_reference.output
    output:
        expand(["{path}/{build}.primary_assembly.genome.{suffix}.bt2"],
                 path = ref_root + "/bt2",
                 build = config['ref']['build'],
                 suffix = ['1', '2', '3', '4', 'rev.1', 'rev.2'])
    threads: 8
    params:
        idx_root = ref_root + "/bt2",
        prefix = config['ref']['build'] + ".primary_assembly.genome"
    shell:
        """
        TEMPDIR=$(mktemp -d -t bt2XXXXXXXXXX)
        FA=$TEMPDIR/temp.fa
        gunzip -c {input} > $FA
        bowtie2-build -f $FA --threads {threads} {params.idx_root}/{params.prefix}
        """

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
    
        """

rule get_rs_fragments:
    input: rules.get_reference.output
    output: rs_frags
    params: enzyme = config['ref']['enzyme']
    threads: 1
    shell:
        """
        # Unzip the genome
        TEMPDIR=$(mktemp -d -t digXXXXXXXXXX)
        FA=$TEMPDIR/temp.fa
        gunzip -c {input} > $FA

        # Get the latest version from the HiC-Pro repo
        wget \
          -O "scripts/digest_genome.py" \
          "https://raw.githubusercontent.com/nservant/HiC-Pro/master/bin/utils/digest_genome.py"

        # Run the python script
        python scripts/digest_genome.py \
          -r {params.enzyme} \
          -o ${output} \
          $FA
        """
