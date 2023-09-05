#!/usr/bin/env Rscript

#===========================================================
# R script for computing the number of TF peaks per anchor
# anchors from differential interactions

# Usage:
#   Rscript Overlap_Peaks_Loops.R --peaks TF_peaks_file.txt --up_loops upregulated_loop_file.txt --down_loops downregulated_loop_file.txt --outdir results/ --prefix Diff_loop_TF_peaks

# Author: Daniela Salgado Figueroa
#===========================================================

library(optparse)
library(foreach)
library(GenomicRanges)
library(ggplot2)
options(scipen=999)

option_list = list(
	make_option(c("--peaks"), default=NULL, type='character', help='TF peaks file with "chr \t start \t end" format'),
	make_option(c("--up_loops"), default=NA, type='character', help='Upregulated loops file with "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format'),
	make_option(c("--down_loops"), default=NA, type='character', help='Downregulated loops file with "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format'),
	make_option(c("--outdir"), default=NULL, type='character', help='Output directory'),
	make_option(c("--prefix"), default=NULL, type='character', help='Output file name')
)
opt_parser <- OptionParser(option_list = option_list, description = "")
opt <- parse_args(opt_parser)

peaks_df <- read.delim(opt$peaks)
loops_up <- read.delim(opt$up_loops)
loops_down <- read.delim(opt$down_loops)
outDir <- opt$outdir
prefix <- opt$prefix

##### Data for test ########### 
#workdir <- '/mnt/bioadhoc-temp/Groups/vd-ay/dsfigueroa/projects/IKAROS'
#datadir <- '/mnt/BioAdHoc/Groups/vd-ay/dsfigueroa/projects/IKAROS'

# Read peaks files
#ikaros_peaks <- read.delim(paste0(datadir, "/March2022/Data/MACS2/MACS2_IKAROS_MAY2020_narrow/IKWTpreBPeaksConc_peaks.narrowPeak"), header=FALSE)

## H3K27ac
#k27ac_diff <- read.delim(paste0(datadir, '/March2022/Figure3/Data/PESLoops/H3K27AC_EdgeR_glmQLFTest_AllLoop_FDR_IHW_DiffLoops_10K_AVP_filtered_Labeled_loops.txt'))
#k27ac_diff$id <- paste0(k27ac_diff$chr1,"-", k27ac_diff$start1,"-", k27ac_diff$end1,"-", k27ac_diff$chr2,"-", k27ac_diff$start2,"-", k27ac_diff$end2)
#k27ac_up <- k27ac_diff[k27ac_diff$logFC < 0, ]
#k27ac_down <- k27ac_diff[k27ac_diff$logFC > 0,]

#peaks_df <- ikaros_peaks
#loops_up <- k27ac_up
#loops_down <- k27ac_down
#outDir <- "./"
#prefix <- "test"

################################################ 

#===============
# Overlap peaks and 
# loop anchors
#===============

loops_list <- list(Up=loops_up, Down=loops_down)

df_loops <- foreach(i=1:length(loops_list), .combine="rbind") %do% {
	
	loops <- loops_list[[i]]
	loop_id <- names(loops_list)[i]

	peaks <- peaks_df
	loops$loop_type <- loop_id
	
	# Overlap peaks and loop anchors
	anch1_gr <- GRanges(loops[,1], IRanges(loops[,2], loops[,3]))
	anch2_gr <- GRanges(loops[,4], IRanges(loops[,5], loops[,6]))
	peaks_gr <- GRanges(peaks[,1], IRanges(peaks[,2], peaks[,3]))
	
	ov_anch1 <- as.data.frame(findOverlaps(peaks_gr, anch1_gr, type="within"))
	ov_anch2 <- as.data.frame(findOverlaps(peaks_gr, anch2_gr, type="within"))
	
	# Extract number of peaks at each anchors
	peaks_anch1 <- as.data.frame(table(ov_anch1$subjectHits))
	peaks_anch1$Var1 <- as.numeric(as.character(peaks_anch1$Var1))

	peaks_anch2 <- as.data.frame(table(ov_anch2$subjectHits))
	peaks_anch2$Var1 <- as.numeric(as.character(peaks_anch2$Var1))

	loops$Anch1_peaks <- 0
	loops$Anch1_peaks[peaks_anch1$Var1] <- peaks_anch1$Freq

	loops$Anch2_peaks <- 0
	loops$Anch2_peaks[peaks_anch2$Var1] <- peaks_anch2$Freq
	
	filename <- paste0(outDir, "/", loop_id, "_", prefix, "_Overlap.txt")
	write.table(loops, file = filename, quote = FALSE, sep = "\t", 
	            row.names = FALSE, col.names = TRUE) 

	# Summarize peak counts 
	anchors_1 <- unique(loops[,c("chr1", "start1", "end1", "Anch1_peaks")])
	anchors_2 <- unique(loops[,c("chr2", "start2", "end2", "Anch2_peaks")])
	colnames(anchors_1) <- c("chr", "start", "end", "peaks")
	colnames(anchors_2) <- c("chr", "start", "end", "peaks")
	
	anchors <- unique(rbind(anchors_1, anchors_2))
	
	df <- data.frame(peaks_count=c("0", "1", "2+"),
	anchors=c(sum(anchors$peaks==0), 
	sum(anchors$peaks==1),
	sum(anchors$peaks>=2)))
	df$percentage=(df$anchors/dim(anchors)[1])*100
	df$Loops <- loop_id

	df
}
df_loops$Loops <- factor(df_loops$Loops, levels=c("Down", "Up"))

#===============
# Plot peak number density
#===============

p <- ggplot(data=df_loops, aes(x=peaks_count, y=percentage, fill=Loops, width=.75)) +
	geom_bar(stat="identity", position=position_dodge()) +
	scale_fill_manual(values=c("#80B1D3", "#FB8072")) +
	scale_y_continuous(limits=c(0, 93)) + 
	labs(x = "", y = "Percentage", title="Peaks per loop anchor") + theme_bw()

pdf(paste0(outDir, "/", prefix, "_Overlap_Barplot.pdf"), height=3, width=3)
print(p)
dev.off()


