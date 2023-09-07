#!/usr/bin/env Rscript

#===========================================================
# R script for annotating anchors into promoters (P), 
# enhancer (E) and structural (S).

# Usage:
#	Rscript PES_Annotation.R --loops loop_file.txt --genes genes_file.gtf --ctcf ctcf_peaks.txt --h3K27ac h3K27ac_peaks.txt --h3k4me1 h3k4me1_peaks.txt --slack 0 --outfile loops_PES_Annotated.pdf

# Author: Sourya Bhattacharyya & Daniela Salgado Figueroa
#===========================================================

library(optparse)
library(GenomicRanges)
library(data.table)
library(readr)
options(scipen=10)
options(warn = -1)
options(datatable.fread.datatable=FALSE)

##================
Overlap1D <- function(Inpdata1, Inpdata2, boundary=1, offset=0, uniqov=TRUE, chrWise=FALSE) {

	if (chrWise == FALSE) {
		ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))
		if (uniqov == TRUE) {
			ov_idx_file1 <- unique(ov1[,1])
			ov_idx_file2 <- unique(ov1[,2])		
		} else {
			ov_idx_file1 <- ov1[,1]
			ov_idx_file2 <- ov1[,2]
		}
		nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1)), ov_idx_file1)
		nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2)), ov_idx_file2)

	} else {
		# get the chromosome list from the given peaks
		ChrList <- union(unique(Inpdata1[,1]), unique(Inpdata2[,1]))

		ov_idx_file1 <- c()
		ov_idx_file2 <- c()
		nonov_idx_file1 <- c()
		nonov_idx_file2 <- c()
		
		# process individual chromosomes in a loop
		for (idx in (1:length(ChrList))) {
			CurrChr <- ChrList[idx]
			# peak indices of first data corresponding to the current chromosome
			Inpdata1_CurrChr_IdxSet <- which(Inpdata1[,1] == CurrChr)
			# peak indices of second data corresponding to the current chromosome
			Inpdata2_CurrChr_IdxSet <- which(Inpdata2[,1] == CurrChr)
			# proceed if both data have some entries for the current chromosome
			if ((length(Inpdata1_CurrChr_IdxSet) == 0) | (length(Inpdata2_CurrChr_IdxSet) == 0)) {
	 			next
			}
			# extract the input datasets limited to the current chromosome
			Inpdata1_CurrChr <- Inpdata1[Inpdata1_CurrChr_IdxSet, ]
			Inpdata2_CurrChr <- Inpdata2[Inpdata2_CurrChr_IdxSet, ]
			
			# compute the overlap
			ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1_CurrChr[,1], IRanges(Inpdata1_CurrChr[,2]+boundary-offset, Inpdata1_CurrChr[,3]-boundary+offset)),GRanges(Inpdata2_CurrChr[,1], IRanges(Inpdata2_CurrChr[,2]+boundary-offset, Inpdata2_CurrChr[,3]-boundary+offset))))
			
			# overlapping / non-overlapping indices for both categories
			if (uniqov == TRUE) {
				temp_ov_idx_file1 <- unique(ov1[,1])
				temp_ov_idx_file2 <- unique(ov1[,2])		
			} else {
				temp_ov_idx_file1 <- ov1[,1]
				temp_ov_idx_file2 <- ov1[,2]
			}
			temp_nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1_CurrChr)), temp_ov_idx_file1)
			temp_nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2_CurrChr)), temp_ov_idx_file2)

			# indices with respect to original loop indices
			# (i.e. with respect to Inpdata1 and Inpdata2)
			if (length(temp_ov_idx_file1) > 0) {
				ov_idx_file1_currchr <- Inpdata1_CurrChr_IdxSet[temp_ov_idx_file1]	
				# append in the final list
				ov_idx_file1 <- c(ov_idx_file1, ov_idx_file1_currchr)
			}
			
			if (length(temp_ov_idx_file2) > 0) {
				ov_idx_file2_currchr <- Inpdata2_CurrChr_IdxSet[temp_ov_idx_file2]
				# append in the final list
				ov_idx_file2 <- c(ov_idx_file2, ov_idx_file2_currchr)
			}
			
			if (length(temp_nonov_idx_file1) > 0) {
				nonov_idx_file1_currchr <- Inpdata1_CurrChr_IdxSet[temp_nonov_idx_file1]
				# append in the final list
				nonov_idx_file1 <- c(nonov_idx_file1, nonov_idx_file1_currchr)
			}
			
			if (length(temp_nonov_idx_file2) > 0) {
				nonov_idx_file2_currchr <- Inpdata2_CurrChr_IdxSet[temp_nonov_idx_file2]
				# append in the final list
				nonov_idx_file2 <- c(nonov_idx_file2, nonov_idx_file2_currchr)
			}
		}	# end chromosomes loop
	}

	# return the overlapping and non-overlapping set of indices
	newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2)
	return(newList)

}

Classify_PES_CS4 <- function(LoopData, RefGTFData, CTCF_PeakData, H3K27AC_PeakData, H3K4me1_PeakData) {
	LoopData$Label_Anch1 <- "O"
	LoopData$Gene_Anch1 <- "NA"
	LoopData$PromType_Anch1 <- "NA"
	LoopData$EnhType_Anch1 <- "NA"

	LoopData$Label_Anch2 <- "O"
	LoopData$Gene_Anch2 <- "NA"
	LoopData$PromType_Anch2 <- "NA"
	LoopData$EnhType_Anch2 <- "NA"

	# Annotate promoter anchors 
	print("Annotate promoter anchors")
	OvTSS_Seg1 <- Overlap1D(LoopData[, 1:3], RefGTFData[, 1:3], boundary=0, offset=slack, uniqov=FALSE)
	OvTSS_Seg2 <- Overlap1D(LoopData[, 4:6], RefGTFData[, 1:3], boundary=0, offset=slack, uniqov=FALSE)

	LoopData$Label_Anch1[OvTSS_Seg1$A_AND_B] <- "P"
	LoopData$Label_Anch2[OvTSS_Seg2$A_AND_B] <- "P"

	## Annotate gene of promoter anchors
	#LoopData$Gene_Anch1[OvTSS_Seg1$A_AND_B] <- RefGTFData$GeneName[OvTSS_Seg1$B_AND_A]
	#LoopData$Gene_Anch2[OvTSS_Seg2$A_AND_B] <- RefGTFData$GeneName[OvTSS_Seg2$B_AND_A]

	LoopData$Gene_Anch1 <- ""
	LoopData$Gene_Anch1[unique(OvTSS_Seg1$A_AND_B)] <- unlist(lapply(unique(OvTSS_Seg1$A_AND_B), function(i){
																genes_idx <- OvTSS_Seg1$B_AND_A[which(OvTSS_Seg1$A_AND_B==i)]
																genes <- unique(RefGTFData$GeneName[genes_idx])
																return(paste0(genes, collapse=",")) }))

	LoopData$Gene_Anch2 <- ""
	LoopData$Gene_Anch2[unique(OvTSS_Seg2$A_AND_B)] <- unlist(lapply(unique(OvTSS_Seg2$A_AND_B), function(i){
																genes_idx <- OvTSS_Seg2$B_AND_A[which(OvTSS_Seg2$A_AND_B==i)]
																genes <- unique(RefGTFData$GeneName[genes_idx])
																return(paste0(genes, collapse=",")) }))

	## Annotate promoters with & without CTCF
	OvCTCF_Seg1 <- Overlap1D(LoopData[, 1:3], CTCF_PeakData[, 1:3], boundary=0, offset=slack, uniqov=FALSE)
	OvCTCF_Seg2 <- Overlap1D(LoopData[, 4:6], CTCF_PeakData[, 1:3], boundary=0, offset=slack, uniqov=FALSE)

	TSS_CTCF_Idx_Seg1 <- intersect(OvTSS_Seg1$A_AND_B, OvCTCF_Seg1$A_AND_B)
	TSS_CTCF_Idx_Seg2 <- intersect(OvTSS_Seg2$A_AND_B, OvCTCF_Seg2$A_AND_B)

	LoopData$PromType_Anch1[OvTSS_Seg1$A_AND_B] <- "CTCF_neg"
	LoopData$PromType_Anch2[OvTSS_Seg2$A_AND_B] <- "CTCF_neg"
	LoopData$PromType_Anch1[TSS_CTCF_Idx_Seg1] <- "CTCF_pos"
	LoopData$PromType_Anch2[TSS_CTCF_Idx_Seg2] <- "CTCF_pos"

	# Annotate active enhancer anchors
	print("Annotate active enhancer anchors")
	OvH3K27ACPeak_Seg1 <- Overlap1D(LoopData[, 1:3], H3K27AC_PeakData[, 1:3], boundary=0, offset=0, uniqov=FALSE)
	EnhIdx_Seg1 <- setdiff(OvH3K27ACPeak_Seg1$A_AND_B, OvTSS_Seg1$A_AND_B)

	OvH3K27ACPeak_Seg2 <- Overlap1D(LoopData[, 4:6], H3K27AC_PeakData[, 1:3], boundary=0, offset=0, uniqov=FALSE)
	EnhIdx_Seg2 <- setdiff(OvH3K27ACPeak_Seg2$A_AND_B, OvTSS_Seg2$A_AND_B)

	LoopData$Label_Anch1[EnhIdx_Seg1] <- "E"
	LoopData$Label_Anch2[EnhIdx_Seg2] <- "E"

	## Annotate active enhancer with & without CTCF
	OvCTCF_Seg1 <- Overlap1D(LoopData[, 1:3], CTCF_PeakData[, 1:3], boundary=0, offset=slack, uniqov=FALSE)
	OvCTCF_Seg2 <- Overlap1D(LoopData[, 4:6], CTCF_PeakData[, 1:3], boundary=0, offset=slack, uniqov=FALSE)

	Enh_CTCF_Idx_Seg1 <- intersect(EnhIdx_Seg1, OvCTCF_Seg1$A_AND_B)
	Enh_CTCF_Idx_Seg2 <- intersect(EnhIdx_Seg2, OvCTCF_Seg2$A_AND_B)

	LoopData$EnhType_Anch1[EnhIdx_Seg1] <- "Active_CTCF_neg"
	LoopData$EnhType_Anch2[EnhIdx_Seg2] <- "Active_CTCF_neg"
	LoopData$EnhType_Anch1[Enh_CTCF_Idx_Seg1] <- "Active_CTCF_pos"
	LoopData$EnhType_Anch2[Enh_CTCF_Idx_Seg2] <- "Active_CTCF_pos"

	# Annotate structural anchors
	print("Annotate structural anchors")
	OvCTCFPeak_Seg1 <- Overlap1D(LoopData[, 1:3], CTCF_PeakData[, 1:3], boundary=0, offset=0, uniqov=FALSE)
	StructIdx_Seg1 <- Reduce(setdiff, list(OvCTCFPeak_Seg1$A_AND_B, OvTSS_Seg1$A_AND_B, EnhIdx_Seg1))

	OvCTCFPeak_Seg2 <- Overlap1D(LoopData[, 4:6], CTCF_PeakData[, 1:3], boundary=0, offset=0, uniqov=FALSE)
	StructIdx_Seg2 <- Reduce(setdiff, list(OvCTCFPeak_Seg2$A_AND_B, OvTSS_Seg2$A_AND_B, EnhIdx_Seg2))

	LoopData$Label_Anch1[StructIdx_Seg1] <- "S"
	LoopData$Label_Anch2[StructIdx_Seg2] <- "S"

	# Annotate poised enhancers
	print("Annotate poised enhancers anchors")
	OvH3K4me1Peak_Seg1 <- Overlap1D(LoopData[, 1:3], H3K4me1_PeakData[, 1:3], boundary=0, offset=0, uniqov=FALSE)
	PoisedEnh_Idx_Seg1 <- Reduce(setdiff, list(OvH3K4me1Peak_Seg1$A_AND_B, OvTSS_Seg1$A_AND_B, EnhIdx_Seg1, StructIdx_Seg1))

	OvH3K4me1Peak_Seg2 <- Overlap1D(LoopData[, 4:6], H3K4me1_PeakData[, 1:3], boundary=0, offset=0, uniqov=FALSE)
	PoisedEnh_Idx_Seg2 <- Reduce(setdiff, list(OvH3K4me1Peak_Seg2$A_AND_B, OvTSS_Seg2$A_AND_B, EnhIdx_Seg2, StructIdx_Seg2))

	LoopData$Label_Anch1[PoisedEnh_Idx_Seg1] <- "E"
	LoopData$Label_Anch2[PoisedEnh_Idx_Seg2] <- "E"

	LoopData$EnhType_Anch1[PoisedEnh_Idx_Seg1] <- "Poised"
	LoopData$EnhType_Anch2[PoisedEnh_Idx_Seg2] <- "Poised"

	return(LoopData)

}

##===============================

option_list = list(
	make_option(c("-l", "--loops"), default=NA, type='character', help='Loop file with "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format'),
	make_option(c("-g", "--genes"), default=NA, type='character', help='GTF annotation file with "chr \t TSS start \t TSS end \t" format'),
	make_option(c("-c", "--ctcf"), default=NA, type='character', help='CTCF Peak file with "chr \t start \t end \t" format'),
	make_option(c("-a", "--h3k27ac"), default=NA, type='character', help='H3K27ac Peak file with "chr \t start \t end \t" format'),
	make_option(c("-m", "--h3k4me1"), default=NA, type='character', help='H3K4me1 Peak file with "chr \t start \t end \t" format'),
	make_option(c("-s", "--slack"), default=NA, type='numeric', help='bp of shift allowed for overlap'),
	make_option(c("-o", "--outfile"), default=NA, type='character', help='Name of output file')
)
opt_parser <- OptionParser(option_list = option_list, description = "")
opt <- parse_args(opt_parser)

LoopData <- as.data.frame(read_delim(opt$loops, "\t", col_name=TRUE))
RefGTFData <- fread(opt$genes, header=T)
CTCF_PeakData <- unique(fread(opt$ctcf, header=F)[, 1:3])
H3K27AC_PeakData <- unique(fread(opt$h3k27ac, header=F)[, 1:3])
H3K4me1_PeakData <- unique(fread(opt$h3k4me1, header=F)[, 1:3])
slack <- opt$slack
outFile <- opt$outfile

df <- Classify_PES_CS4(LoopData, RefGTFData, CTCF_PeakData, H3K27AC_PeakData, H3K4me1_PeakData)

fwrite(df, file=outFile, quote = FALSE, sep = "\t")

cat("\n\nNumber of loops:", dim(LoopData)[1])
cat("\nNumber of annotated loops:", dim(df)[1])
cat("\nAnnotated loops saved in", outFile, "\n\n")


