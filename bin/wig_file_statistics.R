#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(reshape2)
library(wig)

#### Get arguments ##
option_list <- list(
  make_option(c("-i", "--input_path"),    type = "character",   default = ".",
              metavar = "character",      help = "Working directory containing wig files to compare"),
  make_option(c("-b", "--both_strands"),  type = "logical",     default = FALSE, action="store_true",
              metavar = "logical",        help = "If true, both strands are put together for correlation"),
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
both_strands  <- opt$both_strands

chrom_length            <- read.csv(chrom_length, header = FALSE, sep = "\t", row.names = 1)
# This line is only an adaption to the import_wig function which cuts chrom names at the first encountered dot
row.names(chrom_length) <- unlist(lapply(strsplit(row.names(chrom_length), "\\."), `[`, 1))
files                   <- list.files(path = input_path, pattern = "\\.wig$")

if(both_strands){
  group_replicates <- function(file_name) {
    sub("_[^_]*$","", file_name)
  }
  
  # Generate tuples based on the file names. Experiments
  combined_strands <- lapply(split(files,sapply(files,group_replicates)),function(group){
    list(forward = group[grep("_forward",group)],
         reverse = group[grep("_reverse",group)])
  })
}

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
  names(filelist) <- files

  # generate a list of chromosome matrices with all values
  for (chr in chromosomes) {
    # if both strands shall have their correlation calculated together they need to be prepared differently
    #   Here the positions of the reverse strand will have the chromosome length added
    #   Then the values of the reverse strand are appended to the forward strand, leading to a list
    #     that contains values of both strands
    if(both_strands){
      # double length for each chromosome since forward and reverse are combined
      # half column length since each sample comes with 2 files in this case
      chromosome_values[[chr]] <- matrix(nrow = chrom_length[chr, ] * 2, ncol = length(files) / 2, data = c(0))
      
      for (i in seq(1, length(combined_strands))) {
        file_forward <- combined_strands[[i]]$forward
        lines_forward       <- filelist[[file_forward]]$chr == chr
        positions_forward   <- filelist[[file_forward]][lines_forward, ]$pos
        values_forward      <- filelist[[file_forward]][lines_forward, ]$val
        
        file_reverse <- combined_strands[[i]]$reverse
        lines_reverse       <- filelist[[file_reverse]]$chr == chr
        positions_reverse   <- filelist[[file_reverse]][lines_reverse, ]$pos
        values_reverse      <- filelist[[file_reverse]][lines_reverse, ]$val
        
        positions_reverse <- sapply(positions_reverse, function(x) x + chrom_length[chr, ])
        
        
        chromosome_values[[chr]][c(positions_forward,positions_reverse), i] <- c(abs(values_forward),abs(values_reverse))
      }
    } else {
      chromosome_values[[chr]] <- matrix(nrow = chrom_length[chr, ], ncol = length(files), data = c(0))
      # loop through our input files and import the values into our matrix
      for (i in seq(1, length(files))) {
        lines       <- filelist[[i]]$chr == chr
        positions   <- filelist[[i]][lines, ]$pos
        values      <- filelist[[i]][lines, ]$val
        chromosome_values[[chr]][positions, i] <- abs(values)
      }
    }
  }

  # build the result list
  for (chr in chromosomes) {
    chr_dat[[chr]] <- list(name = chr, len = as.numeric(chrom_length[chr, ] * 2), val = chromosome_values[[chr]])
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
  print(appended_values)

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
if(both_strands){
  do_cor_analysis(names(combined_strands), import_wig_files(files))
} else {
  do_cor_analysis(files, import_wig_files(files))
}
