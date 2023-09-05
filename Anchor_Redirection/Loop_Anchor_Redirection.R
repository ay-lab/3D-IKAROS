#!/usr/bin/env Rscript

#===========================================================
# R script for finding redirected anchors
# give a set of loops

# Usage:
#	Rscript Loop_Anchor_Redirection.R --up_loops upregulated_loop_file.txt --down_loops downregulated_loop_file.txt --dataset H3K27ac --outdir results/

# Author: Daniela Salgado Figueroa
#===========================================================

library(optparse)
library(InteractionSet)
library(foreach)
library(ggplot2)
library(ggsignif)

option_list = list(
	make_option(c("--up_loops"), default=NA, type='character', help='Upregulated loops file with "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format'),
	make_option(c("--down_loops"), default=NA, type='character', help='Downregulated loops file with "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format'),
	make_option(c("--dataset"), default=NA, type='character', help='Dataset ID for loops'),
	make_option(c("--outdir"), default=NA, type='character', help='Name of output directory')
)
opt_parser <- OptionParser(option_list = option_list, description = "")
opt <- parse_args(opt_parser)

loops_up <- read.delim(opt$up_loops)
loops_down <- read.delim(opt$down_loops)
dataset_id <- opt$dataset
outDir <- opt$outdir

# Add loop direction and dataset name to loops
loops_up$loop_state <- "Up"
loops_down$loop_state <- "Down"

loops_up$dataset <- dataset_id
loops_down$dataset <- dataset_id

# Associate unique ID to loop anchors
loops_up$anch1_id <- paste0(loops_up[,1], "-", loops_up[,2], "-", loops_up[,3])
loops_up$anch2_id <- paste0(loops_up[,4], "-", loops_up[,5], "-", loops_up[,6])
loops_down$anch1_id <- paste0(loops_down[,1], "-", loops_down[,2], "-", loops_down[,3])
loops_down$anch2_id <- paste0(loops_down[,4], "-", loops_down[,5], "-", loops_down[,6])

# Compute loop size
loops_up$loop_size <- loops_up[,5]-loops_up[,2]
loops_down$loop_size <- loops_down[,5]-loops_down[,2]

#===============
# Perform overlap of 
# up and downregulated loop anchors
#===============

interactions <- list(list(Up=loops_up, Down=loops_down))

slack <- 0
down_up_ov <- foreach(i=1:length(interactions), .combine="rbind") %do% {

	dataset <- names(interactions)[i]
	df1 <- interactions[[i]][[1]]
	df2 <- interactions[[i]][[2]]

	# Perform anchors overlap between up and downregulated loop
	df1_anchors1 <- unique(df1[,c("chr1", "start1", "end1", "loop_state", "dataset")])
	df1_anchors2 <- unique(df1[,c("chr2", "start2", "end2", "loop_state", "dataset")])
	colnames(df1_anchors1) <- c("chr", "start", "end", "loop_state", "dataset")
	colnames(df1_anchors2) <- c("chr", "start", "end", "loop_state", "dataset")

	df2_anchors1 <- unique(df2[,c("chr1", "start1", "end1", "loop_state", "dataset")])
	df2_anchors2 <- unique(df2[,c("chr2", "start2", "end2", "loop_state", "dataset")])
	colnames(df2_anchors1) <- c("chr", "start", "end", "loop_state", "dataset")
	colnames(df2_anchors2) <- c("chr", "start", "end", "loop_state", "dataset")

	df1_anchors <- unique(rbind(df1_anchors1, df1_anchors2))
	df2_anchors <- unique(rbind(df2_anchors1, df2_anchors2))

	df1_anchors_gr <- GRanges(seqnames = df1_anchors$chr, ranges = IRanges(df1_anchors$start-slack, end = df1_anchors$end+slack), strand = '*')
	df2_anchors_gr <- GRanges(seqnames = df2_anchors$chr, ranges = IRanges(df2_anchors$start-slack, end = df2_anchors$end+slack), strand = '*')

	ov_anch <- as.data.frame(findOverlaps(df1_anchors_gr, df2_anchors_gr))

	# Select redirected anchors
	df_res <- rbind(df1_anchors[unique(ov_anch$queryHits),],
	df2_anchors[unique(ov_anch$subjectHits),])

	df_res
}

# Select anchors from downregulated loops
down_ov <- down_up_ov[down_up_ov$loop_state == "Down",]
down_ov$anch_id <- paste0(down_ov[,1], "-", down_ov[,2], "-", down_ov[,3])

#===============
# Select loop involving 
# redirected anchors
#===============

loops_down_ov <- loops_down[loops_down$anch1_id %in% down_ov$anch_id | 
							loops_down$anch2_id %in% down_ov$anch_id, ]

loops_up_ov <- loops_up[loops_up$anch1_id %in% down_ov$anch_id | 
						loops_up$anch2_id %in% down_ov$anch_id, ]

#===============
# Plot loop size distribution
#===============

loops_ov <- rbind(loops_down_ov[,c("loop_size", "loop_state", "dataset")],
				loops_up_ov[,c("loop_size", "loop_state", "dataset")])

loops_ov$loop_size_log10 <- log10(loops_ov$loop_size)

p <- ggplot(loops_ov, aes(x=loop_state, y=loop_size_log10, fill=loop_state)) + 
	geom_boxplot() + facet_wrap(.~dataset) +
	geom_signif(comparisons = list(c("Down", "Up")),
	map_signif_level=TRUE, textsize = 7, vjust = 0.5) +
	scale_fill_manual(values=c("#80B1D3", "#FB8072")) +
	labs(x="", y="log10(Distance)") + theme_bw() + theme(text = element_text(size=24))

pdf(paste0(outDir, "/", dataset_id, "_UpDown_Anchor_Overlap_BoxPlot.pdf"), width=6)
print(p)
dev.off()


