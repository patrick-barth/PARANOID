#!/usr/bin/env Rscript
list.of.packages <- c("optparse")
installed_packages <- list.of.packages %in% rownames(installed.packages())
if(any(installed_packages == FALSE)) {
  install.packages(list.of.packages[!installed_packages], repos = 'http://cran.us.r-project.org')
}
invisible(lapply(list.of.packages, library, character.only = TRUE))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input file in tsv format (output of peak-distance.py)", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file name", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input_file <- opt$input
output_file <- opt$output

parsed_input <- read.csv(input_file, 
                         header = TRUE,
                         sep = "\t")


png(output_file)
  
x <- parsed_input[,1]
y <- parsed_input[,2]
  
  title = paste(sub( 'distances_', '', sub( '\\.tsv$', '', basename(input_file) ) ), sep = " " )
  
  plot(x, y,
       type = "l",
       ylim = c(0, max(y)),
       main = title,
       xlab = "Distance", ylab = "#peaks")
  
  invisible(dev.off())


