#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input file in tsv format (output of peak-distance.py)", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file name", metavar = "character"),
  make_option(c("-t","--type"), type="character",default="png",
              help="File type of output plot", metavar="character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input_file <- opt$input
output_file <- opt$output

parsed_input <- read.csv(input_file, 
                         header = TRUE,
                         sep = "\t")

plot <- ggplot(parsed_input,
       aes(x=distance, y=value)) +
  geom_line(stat="identity") +
    xlab("Peak distance (nt)") +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle=0,hjust=0.5))

ggsave(
  paste(opt$output, ".", opt$type, sep = ""),
  plot = plot,
  device = opt$type
) 

plot_full <- ggplot(parsed_input,
       aes(x=distance, y=value)) +
  geom_line(stat="identity") +
    xlab("Peak distance (nt)") +
    ylim(0,NA) +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle=0,hjust=0.5))

ggsave(
  paste(opt$output, "_full.", opt$type, sep = ""),
  plot = plot_full,
  device = opt$type
)