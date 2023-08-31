#!/usr/bin/env Rscript

#===========================================================
# R script for computing the Aggregate Region Analysis (ARA)
# given a set of interactions

# Usage:
#	Rscript ARA_GENOVA.R --regions loop_file.txt --hic_file contact_map.hic --chr_name 1 --resolution 10000 --min_lim 20 --max_lim 60 --outfile ARA_plot.pdf

# Author: Daniela Salgado Figueroa
#===========================================================

library(optparse)
library(GENOVA)
library(ggplot2)
options(scipen=999)
options(GENOVA.colour.palette="whitered")

option_list = list(
	make_option(c("--regions"), default=NULL, type='character', help='Differential loops (chr \t start \t end'),
	make_option(c("--hic_file"), default=NULL, type='character', help='.hic file for Cat1'),
	make_option(c("--chr_name"), default=1, type='numeric', help='1 for "chrX" format or 0 for "X" format'),
	make_option(c("--resolution"), default=NULL, type='numeric', help='.hic file resolution (ex. 10000)'),
	make_option(c("--min_lim"), default=NULL, type='numeric', help='Min value for color range'),
	make_option(c("--max_lim"), default=NULL, type='numeric', help='Max value for color range'),
	make_option(c("--outfile"), default=NULL, type='character', help='Outfile name')
)
opt_parser <- OptionParser(option_list = option_list, description = "")
opt <- parse_args(opt_parser)

regions <- read.delim(opt$regions)
hicfile <- opt$hic_file
chr_name <- opt$chr_name
res <- opt$resolution
min_limit <- opt$min_lim
max_limit <- opt$max_lim
outFile <- opt$outfile

#===============
# Read contact matrix
#===============

contacts <- load_contacts(hicfile, 
					resolution = res,
					balancing = TRUE,
					sample_name = "")

# Edit name of chromosomes in contact matrix to match loop file
if(chr_name==0){
	contacts$IDX$V1 <- paste0("chr", contacts$IDX$V1)
}

#===============
# Run ARA from GENOVA package
#===============

# size_bin = 100L for 500kb window in each side at 10kb resolution
ara <- ARA(contacts, regions, size_bin = 100L)

if(is.null(max_limit)){
	p1 <- visualise(ara) + 
	theme(text=element_text(size=25), plot.title = element_text(hjust = 0.5)) 
}else{
	p1 <- visualise(ara, colour_lim=c(min_limit, max_limit)) + 
		theme(text=element_text(size=25), plot.title = element_text(hjust = 0.5)) 
}

#===============
# Output ARA plot
#===============

pdf(outFile)
print(p1)
dev.off()



