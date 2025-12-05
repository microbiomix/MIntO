#!/usr/bin/env Rscript

# ---
# title: "Mean read length for strobealign"
# author: "Judit Szarvas"
# date: "2025"
# Estimate mean read length of length-trimmed reads
# ---

# Parse command line arguments
library(optparse)
option_list = list(
  make_option(c("--input"),      type="character", default=NULL, help="qc1 samples read length file", metavar="character"),
  make_option(c("--min_len"),       type="double",    default=NULL, help="minimum length set for trimming", metavar="numeric"),
  make_option(c("--out_meanlen"), type="character", default=NULL, help="output file with mean length value", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (any(is.null(c(opt$input, opt$min_len, opt$out_meanlen)))) {
  print_help(opt_parser)
  stop("Missing required arguments\n", call.=FALSE)
}

# Load libraries
library(dplyr)
library(tidyr)

# Load file with 
read_len_df <- read.table(file=opt$input, header=F, col.names = c('n_reads', 'len_reads', 'sample'), stringsAsFactors = F)

trimmed_min_len <- as.numeric(opt$min_len)

# Filter out read lengths removed by trimming
read_len_df_trimmed <- read_len_df %>% filter(len_reads >= trimmed_min_len)

# Calculate mean length over all
calc_mean_len <- read_len_df_trimmed %>% summarise(mean_len=sum(n_reads*len_reads)/sum(n_reads))

# Corresponding strobealign -r setting as per v0.16.1
calc_mean_len <- calc_mean_len %>% mutate(sba_mean_len=case_when(
  mean_len <=  70 ~  50, 
  mean_len <=  90 ~  75, 
  mean_len <= 110 ~ 100, 
  mean_len <= 135 ~ 125, 
  mean_len <= 175 ~ 150, 
  mean_len <= 375 ~ 250, 
  mean_len <= 500 ~ 400 )
  )

cat(calc_mean_len$sba_mean_len, file=opt$out_meanlen, sep="\n")

