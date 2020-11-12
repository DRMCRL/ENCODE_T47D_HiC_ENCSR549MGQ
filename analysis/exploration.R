library(GenomicInteractions)
library(InteractionSet)
library(rtracklayer)
library(tidyverse)
library(vroom)
library(Gviz)
library(glue)
library(plyranges)

# Let's start working with the standard files (40K)
sz <- 40000
bins <- import.bed(
  glue("data/hic/hic_results/matrix/merged/raw/{sz}/merged_{sz}_abs.bed")
)
sq <- Seqinfo(
  seqnames = seqlevels(bins),
  seqlengths = bins %>% range() %>% end(),
  genome = "GRCh37"
)
seqinfo(bins) <- sq
## The FDR values are wrong and we don't need the bias columns
## Skip those and import everything else. Unfortunately, 3 columns (1:3)
## should be integers, but are doubles
cis_int <- vroom(
  glue("output/MaxHiC/merged/{sz}/cis_interactions.txt.gz"),
  col_types = "dddd-d--"
)

## GenomicInteractions objects have issues using subset
## Use a Ginteractions object instead
gi <- GInteractions(
  anchor1 = gr[as.integer(cis_int$bin1ID)],
  anchor2 = gr[as.integer(cis_int$bin2ID)],
  counts = as.integer(cis_int$observed_interactions),
  exp_interactions = cis_int$exp_interactions,
  p = exp(-cis_int$neg_ln_p_val)
)
seqinfo(gi) <- sq
## Add the correct FDR
gi$fdr = p.adjust(gi$p, "BH")
## Add the distance between pairs
gi$distance <- pairdist(gi)
gi$logRatio <- log2(gi$counts / gi$exp_interactions)
## Export the file as is
write_rds(
  x = gi,
  file = here::here(
    glue("output/gi_{sz}.rds")
  ),
  compress = "gz"
)
# Remove the original file & collect garbage
rm(cis_int)
gc()

# Plot the association between pairwise distances & p-values
# There are a lot of points here
mcols(gi) %>%
  as.data.frame() %>%
  sample_n(1e6) %>%
  ggplot(aes(distance, -log10(p))) +
  geom_point() +
  geom_smooth(se = FALSE) +
  scale_x_continuous(
    labels = scales::comma
  ) +
  theme_bw()

# Check the relationship between p-values & obs/exp
mcols(gi) %>%
  as.data.frame() %>%
  sample_n(1e6) %>%
  ggplot(aes(logRatio, -log10(p))) +
  geom_point()

gi <- read_rds(
  here::here(
    glue("output/gi_{sz}.rds")
  )
)
## Now only keep the significant interactions and collect garbage
gi <- subset(gi, fdr < 0.05)
gc()


# We need the bed files with ChIP peaks to really explore & figure the visualisations out now
# EHF is chr11:34642640-34700000
ehf <- subset(genesGR, gene_name %in% c("EHF"))
# Get the peaks which overlap EHF
e2Peaks <- import.bed("data/external/BED/Consensus_T47D_AR_E2_peaks_only.bed", seqinfo = sq)
dhtPeaks <- import.bed("data/external/BED/Consensus_T47D_AR_E2_DHT_peaks.bed", seqinfo = sq)
ehfPeaks <- subsetByOverlaps(dhtPeaks, ehf)
# See if any have significant interactions
gi %>% subsetByOverlaps(ehfPeaks)
# Map these to genes
genesGR %>%
  subsetByOverlaps(
    gi %>%
      subsetByOverlaps(ehfPeaks)
  )

# Maybe we need a track which is the HiC bin and we can overlay the ChIP Peaks and transcripts with those,
# then add interactions
gen <- "hg19"
chr <- "chr11"
ideo <- IdeogramTrack(genome = gen, chromosome = chr)
ax <- GenomeAxisTrack()
peaks <- ehfPeaks %>%
  AnnotationTrack(name = "ChIP Peaks", col = "black")
gr <- genesGR %>%
  subsetByOverlaps(
    subsetByOverlaps(gi, ehfPeaks)
  ) %>%
  range(ignore.strand = TRUE)
bt <- bins %>%
  subsetByOverlaps(gr) %>%
  AnnotationTrack(name = "HiC Bins")
it <- gi %>%
  subsetByOverlaps(ehfPeaks) %>%
  as("GenomicInteractions") %>%
  InteractionTrack(name = glue("{sz/1000}KB Interactions"))
gm <- geneModels %>%
  subsetByOverlaps(gr) %>%
  GeneRegionTrack(
    name = "Genes",
    transcriptAnnotation = "symbol"
  )
tm <-  transModels %>%
  subsetByOverlaps(gr) %>%
  GeneRegionTrack(
    name = "Transcripts",
    transcriptAnnotation = "symbol"
  )
plotTracks(
  list(ideo, ax, it, peaks, gm)
)
# Now we just need to tidy this up
# We could also just change the genes to genes and add a chromHMM track



# Bed file is
# https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/E119_25_imputed12marks_dense.bed.gz


# To add these to the track use dput(chromCols) & paste the output in
chromCols <- setNames(stateMap$itemRgb, stateMap$name)
hmecHMM <- import.bed("data/external/E119_25_imputed12marks_dense.bed.gz", seqinfo = sq) %>%
  sort()
hmecHMM$cell <- "HMEC"
hmmcHMM <- import.bed("data/external/E027_25_imputed12marks_dense.bed.gz", seqinfo = sq) %>%
  sort()
hmmcHMM$cell <- "HMMC"
hmecTrack <- hmecHMM %>%
  subsetByOverlaps(gr) %>%
  AnnotationTrack(
    name = "HMEC",
    id = .$name,
    stacking = "dense",
    feature = .$name,
    group = .$name,
    col = "transparent",
    `1_TssA` = "#FF0000", `10_TxEnh5'` = "#C2E105", `11_TxEnh3'` = "#C2E105",
    `12_TxEnhW` = "#C2E105", `13_EnhA1` = "#FFC34D", `14_EnhA2` = "#FFC34D",
    `15_EnhAF` = "#FFC34D", `16_EnhW1` = "#FFFF00", `17_EnhW2` = "#FFFF00",
    `18_EnhAc` = "#FFFF00", `19_DNase` = "#FFFF66", `2_PromU` = "#FF4500",
    `20_ZNF/Rpts` = "#66CDAA", `21_Het` = "#8A91D0", `22_PromP` = "#E6B8B7",
    `23_PromBiv` = "#7030A0", `24_ReprPC` = "#808080", `25_Quies` = "#FFFFFF",
    `3_PromD1` = "#FF4500", `4_PromD2` = "#FF4500", `5_Tx5'` = "#008000",
    `6_Tx` = "#008000", `7_Tx3'` = "#008000", `8_TxWk` = "#009600",
    `9_TxReg` = "#C2E105"
    )
hmmcTrack <- hmmcHMM %>%
  subsetByOverlaps(gr) %>%
  AnnotationTrack(
    name = "HMMC",
    id = .$name,
    stacking = "dense",
    feature = .$name,
    group = .$name,
    col = "transparent",
    `1_TssA` = "#FF0000", `10_TxEnh5'` = "#C2E105", `11_TxEnh3'` = "#C2E105",
    `12_TxEnhW` = "#C2E105", `13_EnhA1` = "#FFC34D", `14_EnhA2` = "#FFC34D",
    `15_EnhAF` = "#FFC34D", `16_EnhW1` = "#FFFF00", `17_EnhW2` = "#FFFF00",
    `18_EnhAc` = "#FFFF00", `19_DNase` = "#FFFF66", `2_PromU` = "#FF4500",
    `20_ZNF/Rpts` = "#66CDAA", `21_Het` = "#8A91D0", `22_PromP` = "#E6B8B7",
    `23_PromBiv` = "#7030A0", `24_ReprPC` = "#808080", `25_Quies` = "#FFFFFF",
    `3_PromD1` = "#FF4500", `4_PromD2` = "#FF4500", `5_Tx5'` = "#008000",
    `6_Tx` = "#008000", `7_Tx3'` = "#008000", `8_TxWk` = "#009600",
    `9_TxReg` = "#C2E105"
  )

plotTracks(
  list(ideo, ax, it, hmecTrack, hmmcTrack, peaks, gm)
  # from = start(gr),
  # to = end(gr)
)


