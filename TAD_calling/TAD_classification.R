#!/usr/bin/env Rscript

#===========================================================
# R script for classifying TADs into common, splitted, merged and shifted

# Usage:
#	Rscript TAD_classification.R --tads1 TADs_condition1.txt --tads2 TADs_condition2.txt --cond1 Condition1 --cond2 Condition2 --outdir results/

# Author: Daniela Salgado Figueroa
#===========================================================

library(optparse)
library(InteractionSet)
library(plyr)
options(scipen=999)

##================
Extract_Merge_Split_TADs <- function(tad_set1, tad_set2, cond1, cond2, slack=50000) {

	tad_set1_gr <- GRanges(seqnames=tad_set1[,1], ranges=IRanges(tad_set1[,2], end=tad_set1[,3]))
	tad_set2_gr <- GRanges(seqnames=tad_set2[,1], ranges=IRanges(tad_set2[,2], end=tad_set2[,3]))

	# Classify TADs
	tad_set1$class1 <- ""
	tad_set1[,paste0("tad_", cond2, "_id")] <- ""
	tad_set2$class1 <- ""

	for(i in 1:length(tad_set1[,1])) {
		#cat("\ni=", i)
		tad <- tad_set1[i,]
		tad_gr <- GRanges(seqnames=tad[,1], ranges=IRanges(tad[,2]-slack, end=tad[,3]+slack))

		# Find TADs from set 2 that are within TAD i
		ov <- as.data.frame(findOverlaps(query=tad_set2_gr, subject=tad_gr, type="within"))
		tad_set2_ov <- tad_set2[unique(ov$queryHits),]
		tad_set2_ov <- tad_set2_ov[order(tad_set2_ov[,2], decreasing=FALSE),]

		if( dim(tad_set2_ov)[1] > 0){
			tad_set2_ov$consecutive_flag <- ""

			# Test if TADs from set 2 are consecutive
			for(j in  1:(dim(tad_set2_ov)[1]-1)){
				tad_end <- tad_set2_ov[j, 3]
				idx_consecutive <- which(tad_set2_ov[,2]==tad_end)

				if(length(idx_consecutive) > 0){
					tad_set2_ov$consecutive_flag[j] <- "consecutive"
					tad_set2_ov$consecutive_flag[idx_consecutive] <- "consecutive"
				}
			}

			# If TADs from set 2 are consecutive, then TAD i is classifed as splitted and
			#  TADs from set 2 are classified as merged 
			if(any(grepl("consecutive", tad_set2_ov$consecutive_flag))){
				
				# Classify TAD i from set 1
				tad_set1$class1[i] <- paste0("splitted_in_", cond2)
				
				# Classify TAD i from set 2
				tad_set2_ov_id <- tad_set2_ov$id[tad_set2_ov$consecutive_flag=="consecutive"]
				tad_set1[i, paste0("tad_", cond2, "_id")] <- paste0(tad_set2_ov_id, collapse=",")
				
				tad_set2$class1[tad_set2$id %in% tad_set2_ov_id] <-  paste0("merged_in_", cond1)

			}
		}
	}

	return(list(tad_set1[tad_set1$class==paste0("splitted_in_", cond2),],
				tad_set2[tad_set2$class==paste0("merged_in_", cond1),]))

}

Classify_Shifted_TADs <- function(tad_set1, tad_set2, cond1, cond2) {

	tad_set1_gr <- GRanges(seqnames=tad_set1[,1], ranges=IRanges(tad_set1[,2], end=tad_set1[,3]))
	tad_set2_gr <- GRanges(seqnames=tad_set2[,1], ranges=IRanges(tad_set2[,2], end=tad_set2[,3]))

	# Classify TADs
	tad_set1$class1 <- "shifted_both_sides"
	tad_set1$class2 <- "shifted_both_sides"
	tad_set1[,paste0("tad_", cond2, "_id")] <- ""

	for(i in 1:length(tad_set1[,1])) {
		#cat("\ni=", i)
		
		tad <- tad_set1[i,]
		tad_gr <- GRanges(seqnames=tad[,1], ranges=IRanges(tad[,2]+1, end=tad[,3]-1))

		# Find TADs from set 2 that overlap TAD i 
		ov <- as.data.frame(findOverlaps(tad_set2_gr, tad_gr))
		tad_set2_ov <- tad_set2[unique(ov$queryHits),]
		tad_set2_ov <- tad_set2_ov[order(tad_set2_ov[,2], decreasing=FALSE),]

		# If TADs from set 1 and set 2 share a boundary, then classify TADs if they result in
		# smaller or bigger TADs in set 2 and which boundary is shifted (3' or 5')
		if(tad[,2] %in% tad_set2_ov[,2]) { # TADs that share start coordinate (5' boundary)
			tad_set2_idx <- which(tad_set2_ov[,2]==tad[,2])

			if(tad[,3] > tad_set2_ov[tad_set2_idx, 3]){
				tad_set1$class1[i] <- paste0("shifted_smaller_", cond2)
			}else{
				tad_set1$class1[i] <- paste0("shifted_bigger_", cond2)
			}
			tad_set1$class2[i] <- "shifted_end_boundary"
		}

		if(tad[,3] %in% tad_set2_ov[,3]) {  # TADs that share end coordinate (3' boundary)
			tad_set2_idx <- which(tad_set2_ov[,3]==tad[,3])

			if(tad[,2] < tad_set2_ov[tad_set2_idx, 2]){
				tad_set1$class1[i] <- paste0("shifted_smaller_", cond2)
			}else{
				tad_set1$class1[i] <- paste0("shifted_bigger_", cond2)
			}
			tad_set1$class2[i] <- "shifted_start_boundary"
		}

		# Add ID from TADs in set 2
		tad_set1[i, paste0("tad_", cond2, "_id")] <- paste0(tad_set2_ov$id, collapse=",")
	}

	return(tad_set1)
}

##================

option_list = list(
	make_option(c("--tads1"), default=NA, type='character', help='File with TAD coordinates from condition 1 with "chr \t start \t end" format'),
	make_option(c("--tads2"), default=NA, type='character', help='File with TAD coordinates from condition 2 with "chr \t start \t end" format'),
	make_option(c("--cond1"), default=NA, type='character', help='Name of condition 1'),
	make_option(c("--cond2"), default=NA, type='character', help='Name of condition 2'),
	make_option(c("--outdir"), default=NA, type='character', help='Output directory')
)
opt_parser <- OptionParser(option_list = option_list, description = "")
opt <- parse_args(opt_parser)

tad_cond1 <- read.delim(opt$tads1, header=FALSE)
tad_cond2 <- read.delim(opt$tads2, header=FALSE)
cond1_name <- opt$cond1
cond2_name <- opt$cond2
outDir <- opt$outdir

# Remove TAD calls from chrY
tad_cond1 <- tad_cond1[-grep("Y", tad_cond1[,1]),]
tad_cond2 <- tad_cond2[-grep("Y", tad_cond2[,1]),]

# Calculate TAD size
tad_cond1$size <- tad_cond1[,3]-tad_cond1[,2]
tad_cond2$size <- tad_cond2[,3]-tad_cond2[,2]

#===============
# Extract common TADs: TADs that have
# the same coordinates across the two conditions
#===============
print("Extract common TADs")

tad_cond1$id <- paste0(tad_cond1[,1], "-", tad_cond1[,2], "-", tad_cond1[,3])
tad_cond2$id <- paste0(tad_cond2[,1], "-", tad_cond2[,2], "-", tad_cond2[,3])

ovTAD_cond1 <- which(tad_cond1$id %in% tad_cond2$id)
ovTAD_cond2 <- which(tad_cond2$id %in% tad_cond1$id)

common_cond1 <- tad_cond1[ovTAD_cond1,]
common_cond2 <- tad_cond2[ovTAD_cond2,]

common_cond1$class1 <- "Common"
common_cond2$class1 <- "Common"

# Remove common TADs from TAD sets
tad_cond1 <- tad_cond1[-ovTAD_cond1,]
tad_cond2 <- tad_cond2[-ovTAD_cond2,]

#===============
# Extract splitted and merged TADs.
# Splitted TAD: TADs from one condition the split into 2 or more TADs in the other condition
# Merged TAD: TADs from one condition that get merge into 1 larger TADs in the other condition
#===============
print("Extract splitted and merged TADs")

# Consider a shift of 50kb to classify boundaries overlap
window_size <- 50000

tads_set1_set2 <- Extract_Merge_Split_TADs(tad_set1=tad_cond1, tad_set2=tad_cond2, cond1=cond1_name, cond2=cond2_name, slack=window_size)
tads_set2_set1 <- Extract_Merge_Split_TADs(tad_set1=tad_cond2, tad_set2=tad_cond1, cond1=cond2_name, cond2=cond1_name, slack=window_size)

splitted_cond1 <- tads_set1_set2[[1]]
merged_cond1 <- tads_set2_set1[[2]]

merged_cond2 <- tads_set1_set2[[2]]
splitted_cond2 <- tads_set2_set1[[1]]

# Remove splitted and merged TADs from TAD sets
tad_cond1 <- tad_cond1[!tad_cond1$id %in% c(splitted_cond1$id, merged_cond1$id),]
tad_cond2 <- tad_cond2[!tad_cond2$id %in% c(splitted_cond2$id, merged_cond2$id),]

#===============
# Extract shifted TADs: TADs that share
# only one boundary between the two conditions
#===============
print("Extract shifted TADs")

shifted_cond1 <- Classify_Shifted_TADs(tad_set1=tad_cond1, tad_set2=tad_cond2, cond1=cond1_name, cond2=cond2_name)
shifted_cond2 <- Classify_Shifted_TADs(tad_set1=tad_cond2, tad_set2=tad_cond1, cond1=cond2_name, cond2=cond1_name)

#===============
# Output results
#===============

col_names <- c(names(tad_cond1)[1:3], "size", "id", "class1")

tad_cond1_ann <- Reduce(rbind.fill, list(common_cond1[,col_names], 
						splitted_cond1[,c(col_names, paste0("tad_", cond2_name, "_id"))], 
						merged_cond1[,col_names], 
						shifted_cond1[,c(col_names, "class2", paste0("tad_", cond2_name, "_id"))]))
colnames(tad_cond1_ann)[1:3] <- c("chr", "start", "end")

tad_cond2_ann <- Reduce(rbind.fill, list(common_cond2[,col_names], 
						splitted_cond2[,c(col_names,paste0("tad_", cond1_name, "_id"))], 
						merged_cond2[,col_names], 
						shifted_cond2[,c(col_names, "class2", paste0("tad_", cond1_name, "_id"))]))
colnames(tad_cond2_ann)[1:3] <- c("chr", "start", "end")

write.table(tad_cond1_ann, file=paste0(outDir, "/", cond1_name, "_TADs_annotated.txt"), 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE) 
	
write.table(tad_cond2_ann, file=paste0(outDir, "/", cond2_name, "_TADs_annotated.txt"), 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE) 


