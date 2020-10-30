library(readr)
library(dplyr)
library(rtracklayer)

## These cannot be named but should simply be all matrix files as a tsv,
## Followed by the bed files and the names of the output
args <- commandArgs(TRUE)
# args <- c(
#   ## Matrices
#   "data/hic/hic_results/matrix/ENCLB183QHG/raw/10000/ENCLB183QHG_10000.matrix",
#   "data/hic/hic_results/matrix/ENCLB758KFU/raw/10000/ENCLB758KFU_10000.matrix",
#   ## Bed Files
#   "data/hic/hic_results/matrix/ENCLB183QHG/raw/10000/ENCLB183QHG_10000_abs.bed",
#   "data/hic/hic_results/matrix/ENCLB758KFU/raw/10000/ENCLB758KFU_10000_abs.bed",
#   ## Names of merged outputs
#   "data/hic/hic_results/matrix/merged/raw/10000/merged_10000.matrix",
#   "data/hic/hic_results/matrix/merged/raw/10000/merged_10000_abs.bed"
# )

## Check the two bed files are identical, then just choose one
md5 <- vapply(args[3:4] , tools::md5sum, character(1))
stopifnot(md5[[1]] == md5[[2]])
grl <- args[[3]] %>%
  import.bed %>%
  split(f = seqnames(.))

## Load in both matrix files. This is very RAM intensive as they are big
DF1 <- read_tsv(args[[1]], col_names = c("bin1", "bin2", "count")) %>%
  with(
    DataFrame(
      bin1 = Rle(bin1),
      bin2 = bin2,
      count = count
    )
  )
DF2 <- read_tsv(args[[2]], col_names = c("bin1", "bin2", "count")) %>%
  with(
    DataFrame(
      bin1 = Rle(bin1),
      bin2 = bin2,
      count = count
    )
  )
gc()

## Whilst this could be done using lapply, I can't see any performance gains
## Similarly memory allocation may be a problem, so use a simple for loop
for (i in seq_along(grl)) {

  ## Define the current chromosome & the bin ranges
  x <- grl[[i]]
  message(
    "Merging values for ", as.character(seqnames(x)[1])
  )
  mn <- min(as.numeric(x$name))
  mx <- max(as.numeric(x$name))

  ## Merge the DFs, then collect garbage
  mg <- rbind(
    subset(DF1, bin1 >= mn & bin1 <= mx),
    subset(DF2, bin1 >= mn & bin1 <= mx)
  ) %>%
    sort()
  gc()

  ## Merge the bins & collect garbage
  out <- mg %>%
    as.data.frame() %>%
    group_by(bin1, bin2) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    as("DataFrame")
  gc()

  ## Export the file & collect garbage
  ## Only append for the second chromosome & beyond
  message("Exporting merged data")
  append <- i > 1
  write.table(
    x = out,
    file = args[[5]],
    sep = "\t",
    append = append,
    col.names = FALSE,
    row.names = FALSE
  )

  ## Make sure we have a clean workspace again
  rm(list = c("mg", "out", "x"))
  gc()

}

## Now export the bed file
file.copy(args[[3]], args[[6]], overwrite = TRUE)

