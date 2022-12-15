#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(wig)

option_list <- list(
    make_option(c("-i", "--input_path"),  type = "character",   default = ".",
                metavar = "character",    help = "Tsv file containing strand distribution information"),
    make_option(c("-o", "--output"),      type = "character",   default = "output",
                metavar = "character",    help = "Output file"),
    make_option(c("-t", "--type"),        type = "character",   default = "png",
                metavar = "character",    help = "File type of output plot"),
    make_option(c("-c", "--color"),       type = "character",   default = "#69b3a2",
                metavar = "character",    help = "File type of output plot"),
    make_option(c("-p", "--percentile"),  type = "double",      default = 0.9,
                metavar = "double",       help = "Percentile being used for other analyses")
)
opt_parser  <- OptionParser(option_list=option_list)
opt         <- parse_args(opt_parser)

input_path  <- opt$input_path
output      <- opt$output
type        <- opt$type
color       <- opt$color
percentile  <- opt$percentile / 100

# get all wig files form directory
files <- list.files(path = input_path, pattern = "\\.wig$")

# currently only 2 wig files are supposed to be the input: Error if more or less
if (length(files) > 2) {
  stop("Error: Too many files provided")
}

file_forward <- import_wig(files[1])
file_reverse <- import_wig(files[2])

# get a list of all values and transform them into a data.frame counting the occurences of values
all_peak_values <- c(file_forward$val,abs(file_reverse$val))
value_counts_total <- as.data.frame(table(all_peak_values))
colnames(value_counts_total) <- c("value","Freq")
value_counts_total$value <- as.numeric(as.character(value_counts_total$value))

max_peak <- max(value_counts_total$value)
max_freq <- max(value_counts_total$Freq)
quant <- quantile(all_peak_values, probs = percentile)

# generate plot: Freq*10 & y-axis/10 make values of 1 or less visible
plot <- ggplot(value_counts_total, aes(x = value, y = Freq * 10)) +
  geom_col(width = 1, fill = color) +
  scale_x_continuous(breaks = seq(0, max_peak, 20)) +
  scale_y_log10(labels = function(x) x / 10) +
  geom_vline(xintercept = quant,
             color      = "red",
             linetype   = "dashed") +
  geom_text(aes(x = quant,
                y = max_freq),
            label     = paste("Chosen percentile cutoff (", quant, ")", sep = ""),
            color     = "red",
            hjust     = 0,
            size      = 3,
            position  = position_dodge(width = 1)) +
  labs(title = paste("Peak height distribution of sample ", gsub("_forward.wig","",files[1]), sep = ""),
        x = "Peak height", y = "# of occurences") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(),
        plot.title = element_text(hjust = 0.5))

ggsave(
  paste(output, ".", type, sep = ""),
  plot = plot,
  device = type
)
