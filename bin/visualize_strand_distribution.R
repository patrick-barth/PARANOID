#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Tsv file containing strand distribution information", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="Output file", metavar="character"),
  make_option(c("-t","--type"), type="character",default="png",
              help="File type of output plot", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input       <- opt$input 
output  <- opt$output
type        <- opt$type

file_name     <- tail(unlist(strsplit(input,'/')),n=1)
experiment    <- unlist(strsplit(file_name,'\\.'))[1]
data          <- read.csv(input, header = TRUE, sep = "\t")
chromosomes   <- rep(data[,'chromosome'],
                      each=2)
strand        <- rep(c('forward','reverse'),
                      times=length(row.names(data)))
values        <- data[,c('forward','reverse')]
values        <- as.numeric(unlist(as.list(as.data.frame(t(values)))))
combined_data <- data.frame(chromosomes,strand,values)

plot <- ggplot(combined_data, aes(fill=strand, x=chromosomes, y=values)) +
  geom_bar(position="dodge", stat="identity") +
  labs(title = experiment,
    x='Chromosomes', y='Counts') +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  paste(output, ".", opt$type, sep = ""),
  plot = plot,
  device = opt$type
)
