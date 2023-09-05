#!/usr/bin/env Rscript

#===========================================================
# R script for calling sub-compartments using CALDER 

# Usage:
#	Rscript run_CALDER.R --hic_file contact_map.hic --sample_id condition1 --resolution 10000 --outdir results/

# Author: Daniela Salgado Figueroa
#===========================================================

library(optparse)
library(strawr)
library(CALDER)
source("modified_CALDER_functions.R")
options(scipen=999)

option_list = list(
	make_option(c("--hic_file"), default=NULL, type='character', help='.hic file for Cat1'),
	make_option(c("--sample_id"), default=1, type='character', help='ID of sample'),
	make_option(c("--resolution"), default=NULL, type='numeric', help='.hic file resolution (ex. 10000)'),
	make_option(c("--outdir"), default=NULL, type='character', help='Output directory')
)
opt_parser <- OptionParser(option_list = option_list, description = "")
opt <- parse_args(opt_parser)

hicfile <- opt$hic_file
sample <- opt$sample_id
res <- opt$resolution
outDir <- opt$outdir

#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#genes <- genes(txdb)
#saveRDS(genes, 
#	file = "/home/dsfigueroa/R/x86_64-pc-linux-gnu-library/4.0/CALDER/extdata/TxDb.Mmusculus.UCSC.mm10.knownGene.rds")

# DATA for testing ####################################
datadir <- "/mnt/BioAdHoc/Groups/vd-ay/dsfigueroa/projects/IKAROS/datasets/DNAnexus_JUICER_HIC/HIC/NOVASEQ_JUICER_HIC_STDST_mergedRPL"

hicfile <- paste0(datadir, "/HIC_IKDN_7002_852_Merge_allValidPairsMBOI.hic")
resolution <- 100000
sample <- "HIC_IKDN_7002_852_Merge"
outDir <- "./"
#####################################################

# Call subcompartments for each chromosomes in the mouse genome
for(c in c(1:19, "X")){ 
	cat("\nChromosome: chr", c, "\n")

	# Extract contact count matrix
	contact_mat <- straw(norm="KR", fname=hicfile, 
						chr1loc=paste0("chr", c), chr2loc=paste0("chr", c), unit="BP", binsize=resolution, matrix="observed")
		
	write.table(contact_mat, file=paste0(outDir, "/", sample, "_chr", c, "_matrix_tmp.txt"),
		quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

	dir.create(paste0(outDir, "/", sample))
	# Call subcompartments using the modified version of CALDER
	CALDER_main_v2(paste0(outDir, "/", sample, "_chr", c, "_matrix_tmp.txt"), 
				chr=c, bin_size=resolution, out_dir=paste0(outDir, "/", sample), 
				sub_domains=TRUE, save_intermediate_data=TRUE, genome='mm10')
}

