#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input tsv generated by previous pipeline process", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="out.png", 
              help="Output file name", metavar="character"),
  make_option(c("-t","--type"), type="character",default="png",
              help="File type of output plot", metavar="character"),
  make_option(c("-c", "--color"), type="character", default="#69b3a2", 
              help="Color", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input_file <- opt$input
output_file <- opt$output
color <- opt$color

input_data <- read.table( file = input_file,
                          sep = "\t",
                          header = TRUE)

input_data <- input_data[!(input_data$RNA_subtypes=="total"),]
input_data$RNA_subtypes <- factor(input_data$RNA_subtypes, levels = input_data$RNA_subtypes)

plot <- ggplot(input_data,
      aes(x=RNA_subtypes, y=percentage)) +
  geom_bar(stat="identity", fill = color) +
    ylab("%") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1))

ggsave(
  paste(opt$output, ".", opt$type, sep = ""),
  plot = plot,
  device = opt$type
)