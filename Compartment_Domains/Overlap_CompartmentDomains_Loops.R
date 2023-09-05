#!/usr/bin/env Rscript

#===========================================================
# R script for defining compartment domains 
# diven a set of subcompartments 

# Usage:
#	Rscript Overlap_CompartmentDomains_Loops.R --domain_cond1 Condition1_CompartmentDomains.bed --domain_cond2 Condition2_CompartmentDomains.bed --cond1 Condition1 --cond2 Condition2 --up_loops upregulated_loop_file.txt --down_loops downregulated_loop_file.txt --outdir results/

# Author: Daniela Salgado Figueroa
#===========================================================

library(optparse)
library(GenomicRanges)
library(foreach)
library(ggplot2)

option_list = list(
	make_option(c("--domain_cond1"), default=NULL, type='character', help='Compartment domains outfile from ./Define_CompartmentDomains.R (*_CompartmentDomains.bed) in condition 1'),
	make_option(c("--domain_cond2"), default=NULL, type='character', help='Compartment domains outfile from ./Define_CompartmentDomains.R (*_CompartmentDomains.bed) in condition 2'),
	make_option(c("--cond1"), default=NA, type='character', help='Name of condition 1'),
	make_option(c("--cond2"), default=NA, type='character', help='Name of condition 2'),
	make_option(c("--up_loops"), default=NA, type='character', help='Upregulated loops file with "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format'),
	make_option(c("--down_loops"), default=NA, type='character', help='Downregulated loops file with "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format'),
	make_option(c("--resolution"), default=NULL, type='numeric', help='Interactions resolution (ex. 10000)'),
	make_option(c("--outdir"), default=NULL, type='character', help='Output directory')
)
opt_parser <- OptionParser(option_list = option_list, description = "")
opt <- parse_args(opt_parser)

domain_cond1 <- read.delim(opt$domain_cond1)
domain_cond2 <- read.delim(opt$domain_cond2)
cond1_name <- opt$cond1
cond2_name <- opt$cond2
loops_up <- read.delim(opt$up_loops)
loops_down <- read.delim(opt$down_loops)
resolution <- opt$resolution
outDir <- opt$outdir

#===============
# Overlap loops with compartment domains 
# given a set of interactions
#===============

interactions <- list(loops_down, loops_up)
inter_names <- c("Down", "Up")
domains <- list(domain_cond1, domain_cond2)
dom_names <- c(cond1_name, cond2_name)

slack <- 0

summ_tab <- foreach(i=1:length(interactions), .combine="rbind") %do% { 
	inter <- interactions[[i]]
	inter_id <- inter_names[i]

	dom <- domains[[i]]
	dom_id <- dom_names[[i]]

	# Find compartment domain where midpoint of each loop anchor is located
	anch1_gr <- GRanges(seqnames = inter[,1], ranges = IRanges(inter[,2]+(resolution/2), end = inter[,2]+(resolution/2)), strand = '*')
	anch2_gr <- GRanges(seqnames = inter[,4], ranges = IRanges(inter[,5]+(resolution/2), end = inter[,5]+(resolution/2)), strand = '*')

	comp_gr <- GRanges(seqnames = dom$chr, ranges = IRanges(dom$start-slack, end = dom$end+slack), strand = '*')
	
	ov_anch1 <- as.data.frame(findOverlaps(anch1_gr, comp_gr, type="within"))
	ov_anch2 <- as.data.frame(findOverlaps(anch2_gr, comp_gr, type="within"))

	ov_anchs_df <- as.data.frame(merge(ov_anch1, ov_anch2, by="queryHits", all=TRUE))
	colnames(ov_anchs_df) <- c("inter", "comp_anch1", "comp_anch2")
	ov_anchs_df$class <- abs(ov_anchs_df$comp_anch2-ov_anchs_df$comp_anch1)
	ov_anchs_df$compartment_anch1 <- dom$comp_annotation[ov_anchs_df$comp_anch1]
	ov_anchs_df$compartment_anch2 <- dom$comp_annotation[ov_anchs_df$comp_anch2]


	# Ignore loops with anchors in 2 compartments are inter or intra-compartmental
	ov_anchs_df <- ov_anchs_df[!is.na(ov_anchs_df$class),]
	
	# Classifly interactions into:
	# 	Intra domain (A or B): "class" column = 0
	# 	Inter domain consecutive (A-A, B-A, B-B): "class" column = 1
	# 	Inter domain longerange (A or B): "class" column > 1

	ov_anchs_df$loop_class <- ""
	ov_anchs_df$loop_class[ov_anchs_df$class==0 & ov_anchs_df$compartment_anch1=="A"] <- "intra_A"
	ov_anchs_df$loop_class[ov_anchs_df$class==0 & ov_anchs_df$compartment_anch1=="B"] <- "intra_B"
	ov_anchs_df$loop_class[ov_anchs_df$class==1 & ((ov_anchs_df$compartment_anch1=="A" & ov_anchs_df$compartment_anch2=="B") | 
						(ov_anchs_df$compartment_anch1=="B" & ov_anchs_df$compartment_anch2=="A"))] <- "inter_consecutive"
	ov_anchs_df$loop_class[ov_anchs_df$class>1 & ov_anchs_df$compartment_anch1=="A" & ov_anchs_df$compartment_anch2=="A"] <- "inter_long_AA"
	ov_anchs_df$loop_class[ov_anchs_df$class>1 & ov_anchs_df$compartment_anch1=="B" & ov_anchs_df$compartment_anch2=="B"] <- "inter_long_BB"
	ov_anchs_df$loop_class[ov_anchs_df$class>1 & ((ov_anchs_df$compartment_anch1=="A" & ov_anchs_df$compartment_anch2=="B") | 
						(ov_anchs_df$compartment_anch1=="B" & ov_anchs_df$compartment_anch2=="A"))]  <- "inter_long_AB"

	# Save classification of interactions
	inter[,paste0(dom_id, "_comp_anch1")] <- ""
	inter[,paste0(dom_id, "_comp_anch2")] <- ""
	inter[,paste0(dom_id, "_loop_class")] <- "unclassified"
	inter[ov_anchs_df$inter, paste0(dom_id, "_comp_anch1")] <- ov_anchs_df$compartment_anch1
	inter[ov_anchs_df$inter, paste0(dom_id, "_comp_anch2")] <- ov_anchs_df$compartment_anch2
	inter[ov_anchs_df$inter, paste0(dom_id, "_loop_class")] <- ov_anchs_df$loop_class

	write.table(inter, 
		file = paste0(outDir, "/", inter_id, "_loops_Intra_InterCompartmental_", dom_id, "_Domains.txt"), 
		append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

	# Summarize interaction classes
	df <- data.frame(Class=c(rep("intra", 2), rep("inter_consecutive", 3), rep("inter_longer", 3)),
		Pattern=c("A","B", rep(c("A-*-A", "B-*-B", "B-*-A"),2)), 
		Value=c(sum(ov_anchs_df$loop_class=="intra_A"), 
			sum(ov_anchs_df$loop_class=="intra_B"),
			sum(ov_anchs_df$class==1 & ov_anchs_df$compartment_anch1=="A" & ov_anchs_df$compartment_anch2=="A"),
			sum(ov_anchs_df$class==1 & ov_anchs_df$compartment_anch1=="B" & ov_anchs_df$compartment_anch2=="B"),
			sum(ov_anchs_df$loop_class=="inter_consecutive"),
			sum(ov_anchs_df$loop_class=="inter_long_AA"),
			sum(ov_anchs_df$loop_class=="inter_long_BB"),
			sum(ov_anchs_df$loop_class=="inter_long_AB") ))
	df$Pattern <- factor(df$Pattern, levels=c("A", "B", "A-*-A", "B-*-A", "B-*-B"))
	df$Percentage <- (df$Value/dim(inter)[1])*100
	df$sample <- inter_id 
	df
}

#===============
# Plot proportion of 
# intra- and inter-compartmental loops
#===============

summ_tab$sample <- factor(summ_tab$sample, levels=c("Down", "Up"))
summ_tab$Class[summ_tab$Class == "inter_consecutive"] <- "inter\nconsecutive"
summ_tab$Class[summ_tab$Class == "inter_longer"] <- "inter\nnon-consecutive"

My_Theme = theme(legend.text=element_text(size=20),
			legend.title=element_text(size=20, face="bold"),
			axis.text=element_text(size=20),
			axis.text.x=element_text(size=15, face="bold"),
			axis.title=element_text(size=20),
			strip.text.x = element_text(size = 15))

p1 <- ggplot(data=summ_tab[summ_tab$Class == "intra",], aes(x=sample, y=Percentage, fill=Pattern)) +
		geom_bar(stat="identity", color="black") + 
		scale_fill_manual(values=c("#66C2A5", "#BEBADA", "#FFFFB3")) +
		labs(x="", y="% of differential loops") +
		facet_wrap(~Class, ncol = 2) + theme_bw() + My_Theme

p2 <- ggplot(data=summ_tab[summ_tab$Class != "intra",], aes(x=sample, y=Percentage, fill=Pattern)) +
		geom_bar(stat="identity", color="black") + 
		scale_fill_manual(values=c("#66C2A5", "#BEBADA", "#FFFFB3")) +
		labs(x="", y="% of differential loops") +
		facet_wrap(~Class, ncol = 2) + theme_bw() + My_Theme 

pdf(paste0(outDir, "/Intra_Inter-Compartment_Diff_Loops_BarPlot.pdf"), height=7, width=6)
print(p1)
print(p2)
dev.off()

