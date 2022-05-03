setwd('/Users/taowang9/Documents/2022Spring/ChIP/Motif_mutation')
library(data.table)
library("tidyverse")
library(universalmotif)

args <- commandArgs(trailingOnly=TRUE)
tf <- args[1]

# read the table for all tfs and memes
downloads <- fread('downloads', header=FALSE, col.names=c('tf', 'model', 'narrowPeak', 'motifs'))
memes <- fread('memes.tsv', col.names=c('meme', 'link'), header=FALSE)
memes[, 'tf'] <- downloads[, tf]
memes <- memes %>% column_to_rownames('tf')

# merge data.table
# data <- merge.data.table(downloads, memes, by='tf', sort=FALSE)

# read in the PPM for the corresponding tf and motif
# Note: a tf could have multiple motifs and each motif has one PPM
# tf <- 'E2F1'
meme_prefix <- 'Motif_PPM/'
meme_files <- unlist(strsplit(memes[tf, 'meme'], split=','))

i = 1
for (meme_file in meme_files) {
  meme <- read_meme(paste0(meme_prefix, meme_file))
  mx <- meme@motif
  
  # filter PPM according to pre-defined thresholds for source and target base respectively.
  threshold_source <- 0.75
  threshold_target <- 0.25
  col_max <- apply(mx, 2, max)
  col_min <- apply(mx, 2, min)
  colnames(mx) <- 1:dim(mx)[2]
  mx <- mx[, (col_max > threshold_source) & (col_min < threshold_target)]
  
  # write out the transformation map 
  # index: the index of the source base, 1-based in R, subtract 1 to get python 0-based index
  # source: the source base, in order to transform into target base, the reference base at index position must match source base.
  # target: the target base for the source to transform into
  base_target <- apply(mx, 2, which.min)
  base_source <- apply(mx, 2, which.max)
  target <- rownames(mx)[base_target]
  source <- rownames(mx)[base_source]
  index <- names(base_target)
  probability_source <- apply(mx, 2, max)
  probability_target <- apply(mx, 2, min)
  
  ifelse(!dir.exists(tf), dir.create(tf), FALSE)
  fwrite(data.table(index, target, source, probability_target, probability_source), paste0(tf, "/", tf, '_', 'map_', i, '.csv'), sep='\t')
  i <- i+1
}
