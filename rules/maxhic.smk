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
        mat = "data/hic/hic_results/matrix/{sample}/raw/{bin}/{sample}_{bin}.matrix"
    output:
        cis = "output/MaxHiC/{sample}/{bin}/cis_interactions.txt",
        trans = "output/MaxHiC/{sample}/{bin}/trans_interactions.txt"
    conda: "../envs/maxhic.yml"
    log: "logs/MaxHiC/{sample}_{bin}_MaxHiC.log"
    threads: config['hicpro']['ncpu']
    shell:
        """
        ## Given the problems with the raw output from HiC-Pro, we should
        ## delete the *ord.bed* file
        HICDIR=$(dirname {input.mat})
        OUTDIR=$(dirname {output.cis})

        if compgen -G "$HICDIR/*ord.bed" > /dev/null; then
          echo -e "Deleting unnecessary symlink"
          rm "$HICDIR/*ord.bed""
        fi

        python scripts/MaxHiC/Main.py \
          -t {threads} \
          $HICDIR \
          $OUTDIR &> {log}
        """

