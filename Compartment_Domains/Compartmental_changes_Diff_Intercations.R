#!/usr/bin/env Rscript

#===========================================================
# R script for computing the change of PC values of compartments overlapping 
# anchors from differential interactions

# Usage:
#	Rscript Compartmental_changes_Diff_Intercations.R --comp intra_compartment.bed --cond1 Condition1 --cond2 Condition2 --up_loops upregulated_loop_file.txt --down_loops downregulated_loop_file.txt --outdir results/

# Author: Daniela Salgado Figueroa
#===========================================================

library(optparse)
library(GenomicRanges)
library(ggplot2)
options(scipen=999)

option_list = list(
	make_option(c("--comp"), default=NULL, type='character', help='dcHiC compartments outfile (intra_compartment.bedGraph)'),
	make_option(c("--cond1"), default=NA, type='character', help='Name of condition 1'),
	make_option(c("--cond2"), default=NA, type='character', help='Name of condition 2'),
	make_option(c("--up_loops"), default=NA, type='character', help='Upregulated loops file with "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format'),
	make_option(c("--down_loops"), default=NA, type='character', help='Downregulated loops file with "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format'),
	make_option(c("--outdir"), default=NULL, type='character', help='Output directory')
)
opt_parser <- OptionParser(option_list = option_list, description = "")
opt <- parse_args(opt_parser)

compartments <- read.delim(opt$comp)
cond1_name <- opt$cond1
cond2_name <- opt$cond2
loops_up <- read.delim(opt$up_loops)
loops_down <- read.delim(opt$down_loops)
outDir <- opt$outdir

# Calculate compartmental change 
compartments$comp_change <- compartments[,cond1_name]-compartments[,cond2_name]
compDiff <- compartments[compartments$padj >= (-log10(0.05)), ]
comp_gr <- GRanges(seqnames = compDiff[,1], ranges = IRanges(compDiff[,2], end = compDiff[,3]), strand = '*')

#===============
# Overlap compartment bins and anchors 
# from upregulated and downregulated loops
#===============

interactions <- list(loops_down, loops_up)
inter_names <- c("Down", "Up")

for(i in 1:length(interactions)){
	inter <- interactions[[i]]
	inter_id <- inter_names[i]

	# Extract unique loop anchors
	anch1 <- inter[,1:3]
	anch2 <- inter[,4:6]
	colnames(anch1) <- c("chr", "start", "end")
	colnames(anch2) <- c("chr", "start", "end")
	anchors <- unique(rbind(anch1, anch2))

	# Overlap anchors and compartment bins
	anchors_gr <- GRanges(seqnames = anchors$chr, ranges = IRanges(anchors$start, end = anchors$end), strand = '*')
	ov_anchors <- as.data.frame(findOverlaps(anchors_gr, comp_gr, type="within"))

	# Identify compartment bins that overlap or not overlap anchors 
	comp_ov_anchs <- compDiff[unique(ov_anchors$subjectHits),]
	comp_noOv_anchs <- compDiff[-unique(ov_anchors$subjectHits),]

	comp_ov_anchs$class <- "Anchors Overlapping"
	comp_noOv_anchs$class <- "No Anchors"

	df <- rbind(comp_ov_anchs, comp_noOv_anchs)

	# Plot distribution of compartmental change
	p <- ggplot(df, aes(x=comp_change, fill=class)) +
			geom_density(alpha=.4) +
			labs(title="Compartments and diff interactions", 
				x=paste0("Compartment score difference\n(", cond1_name, " - ", cond2_name, ")"))+
				ylim(0,1.2) +
			theme_bw() +  theme(legend.text=element_text(size=15),
							legend.title=element_text(size=12, face="bold"),
							legend.position="bottom",
							axis.text.y=element_text(size=25),
							axis.text.x=element_text(size=25, face="bold"),
							axis.title=element_text(size=25))

	pdf(paste0(outDir, "/", inter_id, "_Loops_Compartmental_Change_Distribution.pdf"))
	print(p)
	dev.off()

}


