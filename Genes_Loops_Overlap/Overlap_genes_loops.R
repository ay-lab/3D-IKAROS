#!/usr/bin/env Rscript

#===========================================================
# R script for overlapping genes and loop anchors

# Usage:
#	Rscript Overlap_genes_loops.R --genes genes_file.gtf --dgea_file dgea_results.txt --up_loops upregulated_loop_file.txt --down_loops downregulated_loop_file.txt --nd_loops nd_loop_file.txt --outdir results/

# Author: Daniela Salgado Figueroa
#===========================================================

library(optparse)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(rstatix)
library(ggpubr)
library(RColorBrewer)
options(warn = -1)
options(scipen=999)

option_list = list(
	make_option(c("--genes"), default=NA, type='character', help='GTF annotation file with "chr \t TSS start \t TSS end \t" format'),
	make_option(c("--dgea_file"), default=NA, type='character', help='DGEA output file "gene_name \t adj_pvalue \t log2_fold_change" format'),
	make_option(c("--up_loops"), default=NA, type='character', help='Upregulated loops file with "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format'),
	make_option(c("--down_loops"), default=NA, type='character', help='Downregulated loops file with "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format'),
	make_option(c("--nd_loops"), default=NA, type='character', help='Non differential loops file with "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format'),
	make_option(c("--outdir"), default=NA, type='character', help='Name of output directory')
)
opt_parser <- OptionParser(option_list = option_list, description = "")
opt <- parse_args(opt_parser)

refGTFData <-  as.data.frame(fread(opt$genes))
dgea_res <-  as.data.frame(fread(opt$dgea_file))
loops_up <- as.data.frame(fread(opt$up_loops))
loops_down <- as.data.frame(fread(opt$down_loops))
loops_nondiff <- as.data.frame(fread(opt$nd_loops))
outDir <- opt$outdir

# Classify differential and non-differential expressed genes 
dgea_res$gene_state <- "NonDiff"
dgea_res$gene_state[dgea_res$adj_pvalue < 0.05 & dgea_res$log2_fold_change < -1] <- "Downregulated"
dgea_res$gene_state[dgea_res$adj_pvalue < 0.05 & dgea_res$log2_fold_change > 1] <- "Upregulated"

## Save list of down and upregulated genes
write.table(dgea_res$gene_name[dgea_res$gene_state=="Downregulated"],
	file = paste0(outDir, "/Downregulated_DEG.txt"),
	quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE) 

write.table(dgea_res$gene_name[dgea_res$gene_state=="Upregulated"],
	file = paste0(outDir, "/Upregulated_DEG.txt"),
	quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE) 

# Overlap anchors from differential loops with genes
ov_up_anch1 <- as.data.frame(findOverlaps(GRanges(loops_up[,1], IRanges(loops_up[,2], loops_up[,3])),
									GRanges(refGTFData[,1], IRanges(refGTFData[,2], refGTFData[,3]))))
ov_up_anch2 <- as.data.frame(findOverlaps(GRanges(loops_up[,4], IRanges(loops_up[,5], loops_up[,6])),
									GRanges(refGTFData[,1], IRanges(refGTFData[,2], refGTFData[,3]))))

ov_down_anch1 <- as.data.frame(findOverlaps(GRanges(loops_down[,1], IRanges(loops_down[,2], loops_down[,3])),
									GRanges(refGTFData[,1], IRanges(refGTFData[,2], refGTFData[,3]))))
ov_down_anch2 <- as.data.frame(findOverlaps(GRanges(loops_down[,4], IRanges(loops_down[,5], loops_down[,6])),
									GRanges(refGTFData[,1], IRanges(refGTFData[,2], refGTFData[,3]))))

# Extract genes associated with upregulated and downregulated loops
all_genes_up <- c(refGTFData$gene_name[ov_up_anch1$subjectHits], refGTFData$gene_name[ov_up_anch2$subjectHits])
genes_up <- dgea_res[dgea_res$gene_name %in% all_genes_up,]
genes_up <- unique(genes_up)

all_genes_down <- c(refGTFData$gene_name[ov_down_anch1$subjectHits], refGTFData$gene_name[ov_down_anch2$subjectHits])
genes_down <- dgea_res[dgea_res$gene_name %in% all_genes_down,]
genes_down <- unique(genes_down)

# Ignore genes overlapping both upregulated and downregulated loops
genes_up_down <- intersect(genes_up$gene_name, genes_down$gene_name)

genes_up_uniq <- genes_up[!genes_up$gene_name %in% genes_up_down,]
genes_down_uniq <- genes_down[!genes_down$gene_name %in% genes_up_down,]

# Overlap anchors from non-differential loops with genes
refGTFData_nd_loops <- refGTFData[!refGTFData$gene_name %in% c(genes_up$gene_name, genes_down$gene_name),]

ov_nondiff_anch1 <- as.data.frame(findOverlaps(GRanges(loops_nondiff[,1], IRanges(loops_nondiff[,2], loops_nondiff[,3])),
									GRanges(refGTFData_nd_loops[,1], IRanges(refGTFData_nd_loops[,2], refGTFData_nd_loops[,3]))))
ov_nondiff_anch2 <- as.data.frame(findOverlaps(GRanges(loops_nondiff[,4], IRanges(loops_nondiff[,5], loops_nondiff[,6])),
									GRanges(refGTFData_nd_loops[,1], IRanges(refGTFData_nd_loops[,2], refGTFData_nd_loops[,3]))))

# Extract genes associated with non-differential loops
all_genes_nondiff <- c(refGTFData_nd_loops$gene_name[ov_nondiff_anch1$subjectHits], 
					refGTFData_nd_loops$gene_name[ov_nondiff_anch2$subjectHits])
genes_nondiff <- dgea_res[dgea_res$gene_name %in% all_genes_nondiff,]
genes_nondiff <- unique(genes_nondiff)

cat("\nGenes associated with downregulated loops:\n")
table(genes_down_uniq$gene_state)

cat("\nGenes associated with upregulated loops:\n")
table(genes_up_uniq$gene_state)

cat("\nGenes associated with non-differential loops:\n")
table(genes_nondiff$gene_state)

cat("\n")

## Save list of genes associated to down and upregulated loops
write.table(genes_down_uniq$gene_name,
	file = paste0(outDir, "/Genes_Downregulated_Diff_Loops.txt"),
	quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE) 

write.table(genes_up_uniq$gene_name,
	file = paste0(outDir, "/Genes_Upregulated_Diff_Loops.txt"),
	quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE) 

# Plot log2_fold_change of DEG associated with loops

df_loop_gene <- data.frame(loops=c(rep("Upregulated", dim(genes_up_uniq)[1]), 
						rep("Downregulated", dim(genes_down_uniq)[1]), 
						rep("Non-Differential", dim(genes_nondiff)[1]) ),
				gene_type=c(genes_up_uniq$gene_state, genes_down_uniq$gene_state, genes_nondiff$gene_state),
				log2_fold_change=c(genes_up_uniq$log2_fold_change, genes_down_uniq$log2_fold_change, genes_nondiff$log2_fold_change))
df_loop_gene$log2_fold_change[is.na(df_loop_gene$log2_fold_change)] <- 0
df_loop_gene$loops <- factor(df_loop_gene$loops, levels=c("Downregulated", "Upregulated", "Non-Differential"))

## Remove non-differential genes
df_loop_gene <- df_loop_gene[df_loop_gene$gene_type != "NonDiff",]

## Perform statistical test
st_test_deg <- df_loop_gene %>%
	wilcox_test(log2_fold_change ~ loops, alternative = "two.sided", p.adjust.method="bonferroni") %>%
	rstatix::add_significance()
st_test_deg <- st_test_deg %>% add_xy_position(x = "loops")

## Plot
colors <- c(brewer.pal(n = 11, name = "RdBu")[c(10,2)], brewer.pal(n = 8, name = "Set2")[1:2])
My_Theme = theme(legend.text=element_text(size=20),
			axis.text=element_text(size=20),
			axis.text.x=element_text(size=17, face = "bold"),
			axis.title=element_text(size=17),
			legend.title=element_text(size=14))

p1 <- ggviolin(df_loop_gene, x="loops", y="log2_fold_change", fill="loops",
	palette = colors, alpha = 0.7, add = "boxplot", add.params = list(fill = "white")) +
	labs(x="Interactions", y = "Gene log2(Fold Change)") + 
	theme_classic() + My_Theme + theme(legend.position = "none") + 
	stat_pvalue_manual(st_test_deg, bracket.nudge.y = 1, bracket.size = 0.7, size = 8, tip.length = 0.01)

pdf(paste0(outDir, "/Overlap_DEG_Loops_ViolinPlot.pdf"))
print(p1)
dev.off()

