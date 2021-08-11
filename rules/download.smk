rule get_raw_data:
  output:
    expand(
      ["data/raw/fastq/{sample}/{sample}{ext}{suffix}"],
      sample = samples, 
      ext = [read_ext[0], read_ext[1]],
      suffix = suffix
    )
  threads: 1
  shell:
    """
      DEST="data/raw/fastq"

      # The first sample
      wget \
        -O $DEST/ENCLB758KFU/ENCLB758KFU_R1.fastq.gz \
        https://www.encodeproject.org/files/ENCFF639RSW/@@download/ENCFF639RSW.fastq.gz
      wget \
        -O $DEST/ENCLB758KFU/ENCLB758KFU_R2.fastq.gz \
        https://www.encodeproject.org/files/ENCFF901GID/@@download/ENCFF901GID.fastq.gz

      # The second sample
      wget \
        -O $DEST/ENCLB183QHG/ENCLB183QHG_R1.fastq.gz \
        https://www.encodeproject.org/files/ENCFF887JSU/@@download/ENCFF887JSU.fastq.gz
      wget \
        -O $DEST/ENCLB183QHG/ENCLB183QHG_R2.fastq.gz \
        https://www.encodeproject.org/files/ENCFF308DXR/@@download/ENCFF308DXR.fastq.gz
    """
