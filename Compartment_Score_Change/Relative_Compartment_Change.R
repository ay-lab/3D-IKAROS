#!/usr/bin/env Rscript

#===========================================================
# R script for computing the relative compartment change 
# upon induced IKAROS deletion (FigS5F)

# Author: Daniela Salgado Figueroa
#===========================================================

library(ggplot2)
library(foreach)
options(scipen=999)

#===============
# Read dcHiC output
#===============

# Differential compartments across WT, day 0, day 4, day 12 and day 18 after IKAROS induction, and IKDN
comp_induc_all <- read.delim("../Data/dcHiC_WT_IKDN_INDUC_compartments.bedGraph")

# Differential compartments across day 0 and day 18 after IKAROS induction and IKDN
comp_induc_d18 <- read.delim("../Data/dcHiC_IKDN_INDUC_D0_D18_compartments.bedGraph")

# Annotate compartment switches
comp_induc_d18$annotation <- ""
comp_induc_d18$annotation[comp_induc_d18$D0_IKDN > 0 & comp_induc_d18$D18_IKDN < 0 ] <- "AtoB"
comp_induc_d18$annotation[comp_induc_d18$D0_IKDN < 0 & comp_induc_d18$D18_IKDN > 0 ] <- "BtoA"
comp_induc_d18$annotation[comp_induc_d18$D0_IKDN > 0 & comp_induc_d18$D18_IKDN > 0 ] <- "stableA"
comp_induc_d18$annotation[comp_induc_d18$D0_IKDN < 0 & comp_induc_d18$D18_IKDN < 0 ] <- "stableB"

# Extract differential compartments (day 0 vs day 18) PC values in WT, day 4 and day 12after IKAROS induction, and IKDN
comp_d18_diff <- comp_induc_d18[comp_induc_d18$padj >= (-log10(0.05)), ]
compartments <- merge(comp_d18_diff, comp_induc_all, by=c("chr", "start", "end"), all.x=TRUE)
compartments <- compartments[, c("chr", "start", "end", "IKDN", "WT", "D0_IKDN.y", "D4_IKDN", "D12_IKDN", "D18_IKDN.y", "sample_maha.y", "padj.y", "glosh", "annotation")]
colnames(compartments) <- c("chr", "start", "end", "IKDN", "WT", "D0_IKDN", "D4_IKDN", "D12_IKDN", "D18_IKDN", "sample_maha", "padj", "glosh", "annotation")

#===============
# Compute median compartment change 
# of A to B and B to A switches
#===============

st_samples <- c("WT", "D0", "D4", "D12", "D18", "IKDN")
induc_samples <- c("D0", "D4", "D12", "D18")

# Select A to B and B to A compartment switches
compartments <- compartments[compartments$annotation %in% c("AtoB", "BtoA"),]
compartments[is.na(compartments)] <- 0

# Calculate relative change for each compartment bin
relative_change <- lapply(1:dim(compartments)[1], function(i) {

	comp <- compartments[i,]

	max_change_induc <- comp$D18_IKDN-comp$D0_IKDN

	df <- data.frame(idx=i, sample=induc_samples,
		change_induc=c((comp$D0_IKDN-comp$D0_IKDN)/max_change_induc,
			(comp$D4_IKDN-comp$D0_IKDN)/max_change_induc,
			(comp$D12_IKDN-comp$D0_IKDN)/max_change_induc,
			(comp$D18_IKDN-comp$D0_IKDN)/max_change_induc), comp=comp$annotation)
	df
})
change_df <- do.call(rbind, relative_change)
change_df[is.na(change_df)] <- 0
change_df$sample <- factor(change_df$sample, levels=induc_samples)

# Calculate median of relative change in A to B switches
AB_change_median <- data.frame(sample=induc_samples, comp="AtoB")
AB_change_median$median_oe <- unlist(lapply(AB_change_median$sample, function(x) {
							median(change_df$change_induc[change_df$sample==x & change_df$comp=="AtoB"]) }))
AB_change_median$sample <- factor(AB_change_median$sample, levels=st_samples)

# Calculate median of relative change in B to A switches
BA_change_median <- data.frame(sample=induc_samples, comp="BtoA")
BA_change_median$median_oe <- unlist(lapply(BA_change_median$sample, function(x) {
							median(change_df$change_induc[change_df$sample==x & change_df$comp=="BtoA"]) }))
BA_change_median$sample <- factor(BA_change_median$sample, levels=st_samples)

#===============
# Plot
#===============

median_change_df <- rbind(BA_change_median, AB_change_median)
median_change_df$comp <- factor(median_change_df$comp, levels=c("AtoB", "BtoA"))

p <- ggplot(median_change_df, aes(x=sample, y=median_oe, colour=comp)) +
	geom_point(size=5, position=position_dodge(width=0.6)) +
	labs(x="", y="median O/E") +
	scale_colour_manual(values=c("#80B1D3", "#FB8072")) +
	geom_hline(yintercept=1, linetype='dotted') + 
	theme_bw() + theme(text=element_text(size=20)) 

outDir <- "./"
pdf(paste0(outDir, "/Relative_D0_D18_compartment_change.pdf"), width=10)
print(p)
dev.off()


