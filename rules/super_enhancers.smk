rule make_tag_directories:
  input:
    bams = expand(
      ["data/external/H3K27AC/bam/{{sample}}_{rep}_mapped_reads.bam"],
      rep = ['1', '2', '3']
      )
  output:
    tsv = expand(
      ["data/external/H3K27AC/{{sample}}/chr{chr}.tags.tsv"],
      chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']
    ),
    txt = expand(
      ["data/external/H3K27AC/{{sample}}/tag{f}.txt"],
      f = ['Autocorrelation', 'CountDistribution', 'Info', 'LengthDistribution']
    )
  conda: "../envs/homer.yml"
  log: "logs/makeTagDirectory/{sample}.log"
  threads: 2
  shell:
    """
    DIR = $(dirname {output.tsv[0]})
    makeTagDirectory \
      $DIR \
      {input.bams}
    """

rule make_input_tag_directories:
  input:
    bams = "data/external/H3K27AC/bam/T47D_Pooled_input_mapped_reads.bam"
  output:
    tsv = expand(
      ["data/external/H3K27AC/T47D_Pooled_input/chr{chr}.tags.tsv"],
      chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']
    ),
    txt = expand(
      ["data/external/H3K27AC/T47D_Pooled_input/tag{f}.txt"],
      f = ['Autocorrelation', 'CountDistribution', 'Info', 'LengthDistribution']
    )
  conda: "../envs/homer.yml"
  log: "logs/makeTagDirectory/T47D_Pooled_input.log"
  threads: 2
  shell:
    """
    DIR = $(dirname {output.tsv[0]})
    makeTagDirectory \
      $DIR \
      {input.bams}
    """

rule find_super_enhancers:
  input:
    tagdir = os.path.dirname(rules.make_tag_directories.output.tsv[0]),
    indir = os.path.dirname(rules.make_input_tag_directories.output.tsv[0])
  output:
    enh = "data/external/H3K27AC/{sample}/enhancers.txt",
    super_enh = "data/external/H3K27AC/{sample}/superEnhancers.txt"
  conda:  "../envs/homer.yml"
  log: "logs/findPeaks/{sample}.log"
  threads: 2
  shell:
    """
    findPeaks \
      {input.tagdir} \
      -i {input.indir} \
      -style super \
      -typical {output.enh} \
      -o {output.super_enh}
    """
