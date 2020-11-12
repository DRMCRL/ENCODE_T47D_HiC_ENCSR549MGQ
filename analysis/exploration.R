library(GenomicInteractions)
library(InteractionSet)
library(rtracklayer)
library(tidyverse)
library(vroom)
library(Gviz)
library(glue)

# Let's start working with the standard files (40K)
sz <- 40000
gr <- import.bed(
  glue("data/hic/hic_results/matrix/merged/raw/{sz}/merged_{sz}_abs.bed")
)
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

## Now only keep the significant interactions and collect garbage
gi <- subset(gi, fdr < 0.05)
gc()

## Load in the genes. This was downloaded using `wget`
gtf_ftp <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh37_mapping/gencode.v33lift37.annotation.gtf.gz"
gtf <- here::here("data/external/gencode.v33lift37.annotation.gtf.gz")
## Load the complete gtf, then subset to form the useable GR objects
allGR <- import.gff(gtf)
mcols(allGR) <- mcols(allGR) %>%
  .[str_ends(colnames(.), "(id|type|name|status)") & !str_starts(colnames(.), "hgnc|remap|ccds")]
genesGR <- subset(allGR, type == "gene")
mcols(genesGR) <- mcols(genesGR) %>%
  .[str_starts(colnames(.), "gene")]


transGR <- import.gff(gtf, feature.type = "transcript")
%>%
  subset(type %in% c("exon", "UTR"))
mcols(transGR) <- mcols(transGR)[str_ends(colnames(mcols(transGR)), "_(id|type|name|status)")]
transProm <- promoters(transGR) %>%
  GenomicRanges::reduce()

# We need the bed files with ChIP peaks to really explore & figure the visualisations out now
# EHF is chr11:34642640-34700000
ehf <- transGR %>%
  subset(gene_name == "EHF")
# Get the peaks which overlap EHF
e2Peaks <- import.bed("data/external/BED/Consensus_T47D_AR_E2_peaks_only.bed")
ehfPeaks <- subsetByOverlaps(e2Peaks, ehf)
# See if any have significant interactions
gi %>%
  subset(sig) %>%
  subsetByOverlaps(ehfPeaks)
# Map these to genes
genesGR %>%
  subsetByOverlaps(
    gi %>%
      subset(sig) %>%
      subsetByOverlaps(ehfPeaks)
  )

# Maybe we need a track which is the HiC bin and we can overlay the ChIP Peaks and transcripts with those,
# then add interactions
gen <- "hg19"
chr <- "chr11"
ideo <- IdeogramTrack(genome = gen, chromosome = chr)
ax <- GenomeAxisTrack()
fullGM <- granges(allGR)
fullGM$feature <- allGR$type
fullGM$gene <- str_extract(allGR$gene_id, "ENSG[0-9]+")
fullGM$exon <- str_extract(allGR$exon_id, "ENSE[0-9]+")
fullGM$transcript <- str_extract(allGR$transcript_id, "ENST[0-9]+")
fullGM$symbol <- allGR$gene_name
g <- genesGR %>%
  subsetByOverlaps(
    gi %>%
      subset(sig) %>%
      subsetByOverlaps(ehfPeaks)
  ) %>%
  mcols() %>%
  .[["gene_name"]]
bins <- gr %>%
  subsetByOverlaps(
    subset(genesGR, gene_name %in% g) %>%
      range(ignore.strand = TRUE)
  ) %>%
  AnnotationTrack(name = "HiC Bins")
it <- gi %>%
  subset(sig) %>%
  subsetByOverlaps(ehfPeaks) %>%
  as("GenomicInteractions") %>%
  InteractionTrack(name = "40KB Interactions")
peaks <- e2Peaks %>%
  subsetByOverlaps(
    subset(genesGR, gene_name %in% g)
  ) %>%
  AnnotationTrack(name = "ChIP Peaks")
gm <- fullGM %>%
  subsetByOverlaps(
    subset(genesGR, gene_name %in% g) %>%
      range(ignore.strand = TRUE)
  ) %>%
  GeneRegionTrack(
    name = "Genes",
    transcriptAnnotation = "symbol"
  )
plotTracks(
  list(ideo, ax, it, bins, peaks, gm)
)
# Now we just need to tidy this up
# We could also just change the genes to genes and add a chromHMM track
