library(optparse)
library(terra)
library(sf)
library(dplyr)
library(stringr)
library(pbmclapply)
source("/home/lsong/brodie/scripts/calc_conn.R")

option_list <- list(
  make_option(c("-t", "--taxon"),
              action = "store", type = 'character',
              help = "The taxon group in [bird, mammal] to process. [default %default]."),
  make_option(c("-d", "--dists"),
              action = "store", type = 'character',
              help = "The distances to use.")
  )
opt <- parse_args(OptionParser(option_list = option_list))

# Directories and paths
taxon <- opt$taxon
med_dists <- opt$dists
med_dists <- as.integer(strsplit(med_dists, ",")[[1]])
src_dir <- "/home/lsong/brodie/data"
dst_dir <- "/home/lsong/brodie/data"
calc_conn(taxon = taxon, med_dists = med_dists, 
          src_dir = src_dir, dst_dir = dst_dir,
          tosave = FALSE)