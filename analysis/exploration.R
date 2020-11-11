library(GenomicInteractions)
library(InteractionSet)
library(rtracklayer)
library(tidyverse)
library(vroom)
library(Gviz)

# Let's start working with the standard files (40K)
gr <- import.bed("data/hic/hic_results/matrix/merged/raw/40000/merged_40000_abs.bed")
cis_int <- vroom("output/MaxHiC/merged/40000/cis_interactions.txt.gz")
  # mutate(
  #   across(contains("bin"), as.integer),
  #   observed_interactions = as.integer(observed_interactions)
  # )

# gi <- with(
#   cis_int,
#   GenomicInteractions(
#     anchor1 = gr[bin1ID],
#     anchor2 = gr[bin2ID],
#     counts = observed_interactions,
#     p = exp(-neg_ln_p_val),
#     fdr = p.adjust(exp(-neg_ln_p_val), "BH")
#   )
# )
# These objects appear to have trouble with subset
# Try a Ginteractions object instead
gi <- with(
  cis_int,
  GInteractions(
    anchor1 = gr[as.integer(bin1ID)],
    anchor2 = gr[as.integer(bin2ID)],
    counts = as.integer(observed_interactions),
    exp_interactions = exp_interactions,
    p = exp(-neg_ln_p_val),
    fdr = p.adjust(exp(-neg_ln_p_val), "BH")
  )
)
## Add the distance between pairs
gi$distance <- pairdist(gi)
gi$logRatio <- log2(gi$counts / gi$exp_interactions)
gi$sig <- Rle(gi$fdr < 0.05)
## Also remove the pairs with only one count
## as these are always non-significant & will save resources
gi <- subset(gi, counts > 1 & logRatio > 0)

# Remove the original file
rm(cis_int)
gc()

## Check the p-values
hist(gi$p, breaks = 100)
## They look pretty odd to me


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

# Load in the genes
gtf_ftp <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh37_mapping/gencode.v33lift37.annotation.gtf.gz"
gtf <- here::here("data/external/gencode.v33lift37.annotation.gtf.gz")
genesGR <- import.gff(gtf, feature.type = "gene")
transGR <- import.gff(gtf, feature.type = "transcript")
allGR <- import.gff(gtf) %>%
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
