rule make_tag_directories:
  input:
    bams = expand(
      ["data/external/H3K27AC/bam/{{sample}}_{rep}_mapped_reads.bam"],
      rep = ['1', '2', '3']
      )
  output:
    tsv = temp(expand(
      ["data/external/H3K27AC/{{sample}}/chr{chr}.tags.tsv"],
      chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']
    )),
    txt = expand(
      ["data/external/H3K27AC/{{sample}}/tag{f}.txt"],
      f = ['Autocorrelation', 'CountDistribution', 'Info', 'LengthDistribution']
    )
  conda: "../envs/homer.yml"
  log: "logs/makeTagDirectory/{sample}.log"
  threads: 2
  shell:
    """
    DIR=$(dirname {output.tsv[0]})
    makeTagDirectory \
      $DIR \
      {input.bams} &> {log}
    """

rule make_input_tag_directories:
  input:
    bams = "data/external/H3K27AC/bam/T47D_Pooled_input_mapped_reads.bam"
  output:
    tsv = temp(expand(
      ["data/external/H3K27AC/T47D_Pooled_input/chr{chr}.tags.tsv"],
      chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']
    )),
    txt = expand(
      ["data/external/H3K27AC/T47D_Pooled_input/tag{f}.txt"],
      f = ['Autocorrelation', 'CountDistribution', 'Info', 'LengthDistribution']
    )
  conda: "../envs/homer.yml"
  log: "logs/makeTagDirectory/T47D_Pooled_input.log"
  threads: 2
  shell:
    """
    DIR=$(dirname {output.tsv[0]})
    makeTagDirectory \
      $DIR \
      {input.bams} &> {log}
    """

rule find_super_enhancers:
  input:
    tags = rules.make_tag_directories.output.tsv[0],
    input = rules.make_input_tag_directories.output.tsv[0]
  output:
    enh = "output/HOMER/{sample}/enhancers.tsv",
    super_enh = "output/HOMER/{sample}/superEnhancers.tsv"
  conda:  "../envs/homer.yml"
  log: "logs/findPeaks/{sample}.log"
  threads: 2
  shell:
    """
    TAGS=$(dirname {input.tags})
    INPUT=$(dirname {input.input})
    findPeaks \
      $TAGS \
      -i $INPUT \
      -style super \
      -typical {output.enh} \
      -o {output.super_enh} &> {log}
    """
