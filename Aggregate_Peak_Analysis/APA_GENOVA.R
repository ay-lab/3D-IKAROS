#!/usr/bin/env Rscript

#===========================================================
# R script for computing the Aggregate Peak Analysis (APA)
# given a set of interactions

# Usage:
#	Rscript APA_GENOVA.R --loops loop_file.txt --hic_file contact_map.hic --chr_name 1 --resolution 10000 --min_lim 20 --max_lim 60 --outfile APA_plot.pdf

# Author: Daniela Salgado Figueroa
#===========================================================

library(optparse)
library(GENOVA)
library(ggplot2)
options(scipen=999)
options(GENOVA.colour.palette="whitered")

#================================================

# APA score: the value of the center pixel divided by the mean of pixels 15–30 kb downstream of the upstream loci and 15–30 kb upstream (http://boylelab.org/pubs/Bioinformatics_2015_Phanstiel.pdf)
APA_score <- function(APA_Matrix, binsize){
	
	# dimension of output matrix
	nelem <- ncol(APA_Matrix)

	#the middle position of the matrix, corresponding to the interaction region
	Central_row <- as.integer((nelem - 1) / 2) + 1
	Central_col <- Central_row

	# in the MANGO paper, 15-30 Kb (or in general the low and high thresholds provided in this command line option)
	# in positive X and Y directions are used for computing the APA score
	low_pixel_offset <- as.integer(15000 / binsize)
	high_pixel_offset <- as.integer(30000 / binsize)

	APA_score_MANGO_paper <- (APA_Matrix[Central_col, Central_row] * 1.0) / mean(APA_Matrix[(Central_col + low_pixel_offset):(Central_col + high_pixel_offset), (Central_row - low_pixel_offset):(Central_row - high_pixel_offset)])

	return(round(APA_score_MANGO_paper, 2))
}

# APA ratio: compute the ratio of central bin to the remaining matrix
APA_ratio <- function(APA_Matrix){
	
	# dimension of output matrix
	nelem <- ncol(APA_Matrix)

	#the middle position of the matrix, corresponding to the interaction region
	Central_row <- as.integer((nelem - 1) / 2) + 1
	Central_col <- Central_row
	
	APA_Ratio_Central_Rest <- (APA_Matrix[Central_col, Central_row] * 1.0) / mean(c(mean(APA_Matrix[1:(Central_col - 1), 1:(Central_row - 1)]), mean(APA_Matrix[(Central_col + 1):nelem, 1:(Central_row - 1)]), mean(APA_Matrix[1:(Central_col - 1), (Central_row + 1):nelem]), mean(APA_Matrix[(Central_col + 1):nelem, (Central_row + 1):nelem])))

	return(round(APA_Ratio_Central_Rest, 2))
}

#================================================

option_list = list(
	make_option(c("--loops"), default=NULL, type='character', help='File with loops coordinates (chr1 \t start1 \t end1 \t chr2 \t start2 \t end2'),
	make_option(c("--hic_file"), default=NULL, type='character', help='.hic file for Cat1'),
	make_option(c("--chr_name"), default=1, type='numeric', help='1 for "chrX" format or 0 for "X" format'),
	make_option(c("--resolution"), default=NULL, type='numeric', help='.hic file resolution (ex. 10000)'),
	make_option(c("--min_lim"), default=NULL, type='numeric', help='Min value for color range'),
	make_option(c("--max_lim"), default=NULL, type='numeric', help='Max value for color range'),
	make_option(c("--outfile"), default=NULL, type='character', help='Outfile name')
)
opt_parser <- OptionParser(option_list = option_list, description = "")
opt <- parse_args(opt_parser)

loops <- read.delim(opt$loops)
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

# Plot only loops with a distance > 150000
loops$size <- loops[,5]-loops[,2]
loops_plot <- loops[loops$size > 150000,]

# Edit name of chromosomes in contact matrix to match loop file
if(chr_name==0){
	contacts$IDX$V1 <- paste0("chr", contacts$IDX$V1)
}

#===============
# Run APA from GENOVA package
#===============

# size_bin = 11 for 50kb window in each side at 10kb resolution
apa <- APA(contacts, bedpe = loops_plot, size_bin = 11)
apa_score <- APA_score(apply(apa$signal, 2, c), res)
apa_ratio <- APA_ratio(apply(apa$signal, 2, c))

if(is.null(max_limit)){
	p1 <- visualise(apa) + 
		ggtitle(paste0("\nAPA score = ", apa_score, "\nAPA ratio = ", apa_ratio)) + 
		theme(text=element_text(size=25), plot.title = element_text(hjust = 0.5)) 
}else{
	p1 <- visualise(apa, colour_lim=c(min_limit, max_limit)) + 
		ggtitle(paste0("\nAPA score = ", apa_score, "\nAPA ratio = ", apa_ratio)) + 
		theme(text=element_text(size=25), plot.title = element_text(hjust = 0.5)) 
}

#===============
# Output APA plot
#===============

pdf(outFile)
print(p1)
dev.off()


