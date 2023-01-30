#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(reshape2)
library(wig)

#### Get arguments ##
option_list <- list(
  make_option(c("-i", "--input_path"),    type = "character",   default = ".",
              metavar = "character",      help = "Working directory containing wig files to compare"),
  make_option(c("-l", "--chrom_length"),  type = "character",   default = ".",
              metavar = "character",      help = "File containing lengths of chromosomes"),
  make_option(c("-o", "--output"),        type = "character",   default = "output",
              metavar = "character",      help = "Output base name"),
  make_option(c("-t", "--type"),          type = "character",   default = "png",
              metavar = "character",      help = "File type of output plot")
)
opt_parser  <- OptionParser(option_list = option_list)
opt         <- parse_args(opt_parser)

input_path    <- opt$input_path
chrom_length  <- opt$chrom_length
output        <- opt$output
type          <- opt$type

chrom_length            <- read.csv(chrom_length, header = FALSE, sep = "\t", row.names = 1)
# This line is only an adaption to the import_wig function which cuts chrom names at the first encountered dot
row.names(chrom_length) <- unlist(lapply(strsplit(row.names(chrom_length), "\\."), `[`, 1))
files                   <- list.files(path = input_path, pattern = "\\.wig$")

import_wig_files <- function(files = vector())
{
  chromosomes       <- row.names(chrom_length)
  chromosome_values <- list()
  chr_dat           <- list()
  filelist          <- list()

  for (i in seq(1, length(files))) {
    full_path <- paste(input_path, "/", files[i], sep = "")
    filelist[[i]] <- import_wig(full_path)
  }

  # generate a list of chromosome matrices with all values
  for (chr in chromosomes) {
    chromosome_values[[chr]] <- matrix(nrow = chrom_length[chr, ], ncol = length(files), data = c(0))
    # loop through our input files and import the values into our matrix
    for (i in seq(1, length(files))) {
      lines       <- filelist[[i]]$chr == chr
      positions   <- filelist[[i]][lines, ]$pos
      values      <- filelist[[i]][lines, ]$val
      chromosome_values[[chr]][positions, i] <- values
    }
  }

  # build the result list
  for (chr in chromosomes) {
    chr_dat[[chr]] <- list(name = chr, len = chrom_length[chr, ], val = chromosome_values[[chr]])
  }

  return(chr_dat)
}

do_cor_analysis <- function(files = vector(), dat = list()) {
  appended_values <- data.frame(dat[[1]]$val)

  for (chr in seq(2, length(dat))) {
    tmp_df          <- data.frame(dat[[chr]]$val)
    appended_values <- rbind(appended_values, tmp_df)
  }

  appended_values <- appended_values[rowSums(appended_values[]) > 0, ]

  cor_res           <- cor(appended_values)
  colnames(cor_res) <- unlist(lapply(strsplit(basename(files), "_[forward,reverse]+\\."), `[`, 1))
  rownames(cor_res) <- colnames(cor_res)

  csv_output_file <- paste(output, "_correlation.csv", sep = "")
  cat("Writing tabular output to ", csv_output_file, "...\n")
  write.csv(cor_res, file = csv_output_file)

  print(cor_res)
  molten_cor_res <- melt(cor_res)


  # print all the stuff
  plot <- ggplot(molten_cor_res, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
    labs(title = paste("Correlation of samples belonging\nto ", basename(output), sep = "")) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.y = element_blank())

  ggsave(
    paste(output, "_correlation.", type, sep = ""),
    plot = plot,
    device = type
  )
}

do_cor_analysis(files, import_wig_files(files))
