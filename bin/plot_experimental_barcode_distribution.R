#!/usr/bin/env Rscript
#### Package management ##
if( !require("BiocManager", character.only=TRUE )){install.packages("BiocManager");library("BiocManager")}
CRAN_packages <- c("optparse", "ggplot2")
bioconductor_packages <- c()

for(p in CRAN_packages){
  if( !require( p, character.only=TRUE )){ install.packages(p,repos='http://cran.us.r-project.org') }
  library( p, character.only=TRUE )
}
for(p in bioconductor_packages){
  if( !require( p, character.only=TRUE )){ BiocManager::install(p) }
  library( p, character.only=TRUE )
}
#### Get arguments ##
option_list = list(
  make_option(c("-l","--logs"), type="character",default=NULL,
              help="Log file containing information about barcode splitting", metavar="character"),
  make_option(c("-o","--output"), type="character",default="barcode_distribution",
              help="Name of output plot", metavar="character"),
  make_option(c("-t","--type"), type="character",default="png",
              help="File type of output plot", metavar="character"),
  make_option(c("-c", "--color"), type="character", default="#69b3a2", 
              help="Color", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

color <- opt$color

log_file <- read.csv(normalizePath(opt$logs),
                     sep = "\t")
log_file <- log_file[!(log_file$Barcode=="total"),]

plot <- ggplot(log_file,
       aes(x=Barcode, y=Count)) +
  geom_bar(stat="identity", fill = color) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1))

ggsave(
  paste(opt$output, ".", opt$type, sep = ""),
  plot = plot,
  device = opt$type
)
