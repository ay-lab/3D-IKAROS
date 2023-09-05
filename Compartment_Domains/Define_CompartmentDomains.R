#!/usr/bin/env Rscript

#===========================================================
# R script for defining compartment domains 
# diven a set of subcompartments 

# Usage:
#	Rscript Define_CompartmentDomains.R --subcomp subcompartments_cond1.txt --sample_id condition1 --outdir results/

# Author: Daniela Salgado Figueroa
#===========================================================

library(optparse)
library(foreach)
options(scipen=999)

option_list = list(
	make_option(c("--subcomp"), default=NULL, type='character', help='Subcompartment outfile from CALDER (*_sub_compartments.bed)'),
	make_option(c("--sample_id"), default=1, type='character', help='ID of sample'),
	make_option(c("--outdir"), default=NULL, type='character', help='Output directory')
)
opt_parser <- OptionParser(option_list = option_list, description = "")
opt <- parse_args(opt_parser)

subcomp_df <- read.delim(opt$subcomp)
colnames(subcomp_df) <- c("chr", "start", "end", "subcomp", "val", "strand", "start", "end", "color")

sample <- opt$sample_id
outDir <- opt$outdir

#===============
# Extract domains in each chromosome
#===============

domains_sample <- foreach(c=unique(subcomp_df$chr), .combine="rbind") %do% {
	
	subcomp_chr <- subcomp_df[subcomp_df$chr==c,]
	subcomp_chr <- subcomp_chr[order(subcomp_chr$start, decreasing=FALSE),]
	subcomp_count <- dim(subcomp_chr)[1]
	#rownames(subcomp_chr) <- 1:subcomp_count

	domain_id <- 1
	subcomp_chr$domain <- ''
	for(i in 1:subcomp_count){
		
		# Extract compartment type (A or B) of current subcompartment bin
		currbin_comp <- unlist(strsplit(subcomp_chr$subcomp[i], ""))[1]

		# Extract compartment type (A or B) of next subcompartment bin
		# If current bin is the last bin of chromosome, then
		# extract compartment type (A or B) of previous subcompartment bin

		if(i==subcomp_count){ # Assign last and last-1 with last-2 domain
			evalbin_comp <- unlist(strsplit(subcomp_chr$subcomp[i-1], ""))[1]
		}else{
			evalbin_comp <- unlist(strsplit(subcomp_chr$subcomp[i+1], ""))[1]
		}
		
		# Assign same domain_id to consecutive subcompartment bins that are
		# the same compartment type (A or B)
		if(currbin_comp==evalbin_comp){ 
			subcomp_chr$domain[i] <- domain_id 
		}else{
			subcomp_chr$domain[i] <- domain_id 
			domain_id <- domain_id + 1
		}
	}

	# Extract coordinates of compartmental domains
	domain_chr <- foreach(d=unique(subcomp_chr$domain), .combine="rbind") %do%{
		comp_dom <- subcomp_chr[subcomp_chr$domain==d,]
		comp_dom <- comp_dom[order(comp_dom$start, decreasing=FALSE),]

		data.frame(chr=c, start=comp_dom$start[1], end=comp_dom$end[dim(comp_dom)[1]], 
				comp_annotation=unlist(strsplit(comp_dom$subcomp[1], ""))[1] )
	}	
	
	domain_chr
}

#===============
# Output compartment domains
#===============

write.table(domains_sample, 
	file = paste0(outDir, "/", sample, "_CompartmentDomains.bed"),
	quote = FALSE, sep = "\t", row.names = FALSE)





