rule get_maxhic:
    output:
        dir = directory("scripts/MaxHiC")
    threads: 1
    shell:
        """
        wget https://github.com/Rassa-Gvm/MaxHiC/archive/master.zip
        unzip master.zip -d MaxHiC
        mv MaxHiC/MaxHiC-master scripts/MaxHiC
        rm master.zip
        rmdir MaxHiC
        """

rule run_maxhic:
    input:
        maxhic_dir = rules.get_maxhic.output.dir,
        matrix_dir = "data/hic/hic_results/matrix/ENCLB183QHG/raw/40000/",
        bed_file = "data/hic/hic_results/matrix/ENCLB183QHG/raw/40000/ENCLB183QHG_40000_abs.bed",
        interaction_file = "data/hic/hic_results/matrix/ENCLB183QHG/raw/40000/ENCLB183QHG_40000.matrix"
    output:
        test = "output/maxhic.txt"
    conda: "../envs/maxhic.yml"
    shell:
        """
        python scripts/MaxHiC/Main.py -h > {output}
        """

