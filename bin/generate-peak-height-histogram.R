#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)

option_list <- list(
    make_option(c("-i", "--input"),     type="character",   default=NULL,
                metavar="character",    help="Tsv file containing strand distribution information"),
    make_option(c("-o", "--output"),    type="character",   default="output",
                metavar="character",    help="Output file"),
    make_option(c("-t","--type"),       type="character",   default="png",
                metavar="character",    help="File type of output plot"),
    make_option(c("-t","--color"),      type="character",   default="png",
                metavar="character",    help="File type of output plot")
)
opt_parser  <- OptionParser(option_list=option_list)
opt         <- parse_args(opt_parser)

input       <- opt$input 
output      <- opt$output
type        <- opt$type
color       <- opt$color