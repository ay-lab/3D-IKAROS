#!/usr/bin/env Rscript

#===========================================================
# R script for filtering TADs based on differences in contacts between the two conditions

# Usage:
#	Rscript TAD_filtering.R --tads1 condition1_TADs_annotated.txt --tads2 condition2_TADs_annotated.txt --hic_file1 contacts_cond1.hic --hic_file2 contacts_cond2.hic --cond1 Condition1 --cond2 Condition2 --outdir results/

# Author: Daniela Salgado Figueroa
#===========================================================

library(optparse)
library(strawr)
library(foreach)
library(GENOVA)
library(reshape2)
library(ggplot2)
library(tidyr)
options(scipen=999)

##================
Score_Merge_Split_TADs <- function(tad_set1, tad_set2, hic_set1, hic_set2, cond1, cond2) {

	# Select TADs from condition 1 that get splitted in condition 2  
	splitted_tads <- tad_set1[grep("splitted", tad_set1$class1),]

	# Calculate score for merged and splitted TADs
	merged_scored <- foreach(x=1:length(splitted_tads[,1]), .combine="rbind") %do% {
		#print(x)
		tad <- splitted_tads[x,]
		
		# Find TADs from set 2 that are within TAD x
		tads_set2 <- data.frame(id=unlist(strsplit(tad[,paste0("tad_", cond2, "_id")], ",")))
		tads_set2 <- merge(tads_set2, tad_set2, all.x=TRUE)
		tads_set2 <- tads_set2[order(tads_set2$start, decreasing=FALSE),]

		# Find 3' TADs
		last_split_tad <- tads_set2[dim(tads_set2)[1],]

		# Extract score of each upstream Z regions
		z_score <- list()
		run_z_score <- foreach(i=1:(dim(tads_set2)[1]-1)) %do% {
			z_cond1_m <- straw(norm="KR", fname=hic_set1, 
					chr1loc=paste0(tads_set2$chr[i], ":", tads_set2$end[i], ":", last_split_tad$end), 
					chr2loc=paste0(tads_set2$chr[i], ":", tads_set2$start[i], ":", tads_set2$end[i]), 
					unit="BP", binsize=resolution, matrix=matrix_value)

			z_cond2_m <- straw(norm="KR", fname=hic_set2, 
					chr1loc=paste0(tads_set2$chr[i], ":", tads_set2$end[i], ":", last_split_tad$end), 
					chr2loc=paste0(tads_set2$chr[i], ":", tads_set2$start[i], ":", tads_set2$end[i]), 
					unit="BP", binsize=resolution, matrix=matrix_value)

			# Fill all possible values of matrix (give value=0 for misiing values)
			x_coord <- seq(tads_set2$start[i], tads_set2$end[i], resolution)
			y_coord <- seq(tads_set2$end[i], last_split_tad$end, resolution)

			tmp <- data.frame(y=sort(rep(y_coord, length(x_coord))))
			tmp$x <- x_coord

			z_cond1_m_full <- merge(z_cond1_m, tmp, by=c("x", "y"), all=TRUE)
			z_cond1_m_full$counts[is.na(z_cond1_m_full$counts)] <- 0

			z_cond2_m_full <- merge(z_cond2_m, tmp, by=c("x", "y"), all=TRUE)
			z_cond2_m_full$counts[is.na(z_cond2_m_full$counts)] <- 0

			# Calculate median of Z regions
			z_score[[paste0("z_", i)]] <- c(median(z_cond1_m_full$counts), 
										median(z_cond2_m_full$counts))
			names(z_score[[paste0("z_", i)]]) <- c(paste0(cond1, "_median"), paste0(cond2, "_median"))

		}

		# Extract score of TAD region (X region)
		tads_set2[,paste0("TAD_score_", cond1, "_median")] <- NA
		tads_set2[,paste0("TAD_score_", cond2, "_median")] <- NA
		
		tads_set2[,paste0("Z_score_", cond1, "_median")] <- NA
		tads_set2[,paste0("Z_score_", cond2, "_median")] <- NA

		run <- foreach(i=1:(dim(tads_set2)[1])) %do% {
			x_cond1_m <- straw(norm="KR", fname=hic_set1, 
				chr1loc=paste0(tads_set2$chr[i], ":", tads_set2$start[i], ":", tads_set2$end[i]), 
				chr2loc=paste0(tads_set2$chr[i], ":", tads_set2$start[i], ":", tads_set2$end[i]), 
				unit="BP", binsize=resolution, matrix=matrix_value)

			x_cond2_m <- straw(norm="KR", fname=hic_set2, 
				chr1loc=paste0(tads_set2$chr[i], ":", tads_set2$start[i], ":", tads_set2$end[i]), 
				chr2loc=paste0(tads_set2$chr[i], ":", tads_set2$start[i], ":", tads_set2$end[i]), 
				unit="BP", binsize=resolution, matrix=matrix_value)

			# Fill all possible values of matrix
			coord <- seq(tads_set2$start[i], tads_set2$end[i], resolution)
			tmp <- as.data.frame(t(combn(coord, 2)))
			colnames(tmp) <- c("x", "y")
				
			x_cond1_m_full <- merge(x_cond1_m, tmp, all=TRUE)
			x_cond1_m_full$counts[is.na(x_cond1_m_full$counts)] <- 0

			x_cond2_m_full <- merge(x_cond2_m, tmp, all=TRUE)
			x_cond2_m_full$counts[is.na(x_cond2_m_full$counts)] <- 0

			# Calculate median of X regions
			tads_set2[i, paste0("TAD_score_", cond1, "_median")] <- median(x_cond1_m_full$counts)
			tads_set2[i, paste0("TAD_score_", cond2, "_median")] <- median(x_cond2_m_full$counts)

			# Associate median of Z regions
			if(i==1){ 
				tads_set2[i, paste0("Z_score_", cond1, "_median")] <- z_score[[i]][paste0(cond1, "_median")]
				tads_set2[i, paste0("Z_score_", cond2, "_median")] <- z_score[[i]][paste0(cond2, "_median")]
			} else{
				tads_set2[i, paste0("Z_score_", cond1, "_median")] <- z_score[[i-1]][paste0(cond1, "_median")]
				tads_set2[i, paste0("Z_score_", cond2, "_median")] <- z_score[[i-1]][paste0(cond2, "_median")]
			}
		}

		# Calculate score for TAD that get splitted
		splitted_tads[x, paste0("sum_Z_score_", cond1, "_median")] <- sum(unique(tads_set2[paste0("Z_score_", cond1, "_median")]))
		splitted_tads[x, paste0("sum_Z_score_", cond2, "_median")] <- sum(unique(tads_set2[paste0("Z_score_", cond2, "_median")]))
		splitted_tads$score_diff[x] <- splitted_tads[x, paste0("sum_Z_score_", cond1, "_median")]-splitted_tads[x, paste0("sum_Z_score_", cond2, "_median")]

		# Calculate score for TAD that get merged
		tads_set2$score_diff <- tads_set2[,paste0("Z_score_", cond1, "_median")]-tads_set2[,paste0("Z_score_", cond2, "_median")]

		tads_set2
	}
	splitted_tads$score_diff[is.na(splitted_tads$score_diff)] <- 0

	return(list(splitted_tads, merged_scored)) 

}

Score_Shift_TADs <- function(tad_set1, tad_set2, hic_set1, hic_set2, cond1, cond2, slack=50000) {

	# Select shifted TADs that share one boundary with the other condition 
	shifted <- tad_set1[grep("shifted_bigger|shifted_smaller", tad_set1$class1),]

	# Calculate score for shifted TADs
	shifted_scored <- foreach(x=1:dim(shifted)[1], .combine="rbind") %do% {
		#print(x)

		tad_shift_set1 <- shifted[x,]
		tads_set2_id <- data.frame(id = unlist(strsplit(tad_shift_set1[,paste0("tad_", cond2, "_id")], ",")) )
		tads_set2 <- separate(data=tads_set2_id, col=id, into=c("chr", "start", "end"), sep="-")
		
		# Identify shifted TAD from condition 2
		if(tad_shift_set1$start %in% tads_set2$start){
			tad_set2 <- tads_set2[grep(tad_shift_set1$start, tads_set2$start),]
		}else{
			tad_set2 <- tads_set2[grep(tad_shift_set1$end, tads_set2$end),]
		}

		# Define coordinates of differential region to calculate difference in contacts
		if(grepl("bigger", tad_shift_set1$class1)){
			region <- data.frame(chr=tad_set2[1], start=as.numeric(tad_set2[2]), end=as.numeric(tad_set2[3]))
			
			if( grepl("end", tad_shift_set1$class2)) {
				region$midpoint <- region$end
				region$extension <- tad_shift_set1$end
				chr1loc_window <- paste0(region$chr, ":", region$extension, ":", region$midpoint+slack)
			}else{
				region$midpoint <- region$start
				region$extension <- tad_shift_set1$start
				chr1loc_window <- paste0(region$chr, ":", region$midpoint-slack, ":", region$extension)
			}
		}else{
			region <- data.frame(chr=tad_shift_set1[1], start=tad_shift_set1[2], end=tad_shift_set1[3])
			
			if( grepl("end", tad_shift_set1$class2)) {
				region$midpoint <- region$end
				region$extension <- tad_set2$end
				chr1loc_window <- paste0(region$chr, ":", region$extension, ":", region$midpoint+slack)
			}else{
				region$midpoint <- region$start
				region$extension <- tad_set2$start
				chr1loc_window <- paste0(region$chr, ":", region$midpoint-slack, ":", region$extension)
			}
		}

		# Extract contact count of differential regions from both conditions
		m_cond1 <- straw(norm="KR", fname=hic_set1, 
						chr1loc=chr1loc_window,
						chr2loc=paste0(region$chr, ":", region$start, ":", region$end), 
						unit="BP", binsize=resolution, matrix=matrix_value)

		m_cond2 <- straw(norm="KR", fname=hic_set2, 
						chr1loc=chr1loc_window,
						chr2loc=paste0(region$chr, ":", region$start, ":", region$end), 
						unit="BP", binsize=resolution, matrix=matrix_value)

		# Fill all possible values of matrix (give value=0 for misiing values)
		x_coord <- seq(region$start,  region$end, resolution)
		y_coord <- seq(region$midpoint-slack, region$midpoint+slack, resolution)

		tmp <- data.frame(y=sort(rep(y_coord, length(x_coord))))
		tmp$x <- x_coord

		m_cond1_full <- merge(m_cond1, tmp, by=c("x", "y"), all=TRUE)
		m_cond1_full$counts[is.na(m_cond1_full$counts)] <- 0

		m_cond2_full <- merge(m_cond2, tmp, by=c("x", "y"), all=TRUE)
		m_cond2_full$counts[is.na(m_cond2_full$counts)] <- 0

		tad_shift_set1[,paste0("score_X_", cond1)] <- median(m_cond1_full$counts)
		tad_shift_set1[,paste0("score_X_", cond2)] <- median(m_cond2_full$counts)
		
		if(grepl("bigger", tad_shift_set1$class1)){
			tad_shift_set1$score_diff <- tad_shift_set1[,paste0("score_X_", cond2)]-tad_shift_set1[,paste0("score_X_", cond1)]
		}else{
			tad_shift_set1$score_diff <- tad_shift_set1[,paste0("score_X_", cond1)]-tad_shift_set1[,paste0("score_X_", cond2)]
		}

		tad_shift_set1
	}
	shifted_scored$score_diff[is.na(shifted_scored$score_diff)] <- 0

	return(shifted_scored)
}

ContactMapDiff <- function(mat1, mat2, colour_lim){

	# Calculate difference by subtracting contact counts 
	diff <- mat1 - mat2
	lfc <- log2(mat1/mat2)
	df <- reshape2::melt(diff)
	colnames(df) <- c("y", "x", "name", "value")
	df$Difference <- df$value
	#df$Difference[df$value < -15] <- -15

	# Limit color scale
	df$Difference[df$Difference > colour_lim] <- colour_lim
	df$Difference[df$Difference < -colour_lim] <- -colour_lim

	# Plot
	p3 <- ggplot(df, aes(x = x, y = y, fill = Difference)) + 
		geom_tile() +
		scale_fill_gradient2(low = 'forestgreen', high = 'purple', mid = 'white', midpoint = 0, limit=c(-colour_lim, colour_lim)) +
		scale_x_continuous(name="", labels=c("", "5' border", "", "3' border", "")) + 
		scale_y_continuous(name="", labels=c("", "3' border", "", "5' border", "")) +
		theme_minimal() + theme(panel.grid.minor = element_blank())
	
	return(p3)
}

##================

option_list = list(
	make_option(c("--tads1"), default=NA, type='character', help='File with TAD coordinates from condition 1 with "chr \t start \t end" format'),
	make_option(c("--tads2"), default=NA, type='character', help='File with TAD coordinates from condition 2 with "chr \t start \t end" format'),
	make_option(c("--hic_file1"), default=NA, type='character', help='.hic file for condition 1'),
	make_option(c("--hic_file2"), default=NA, type='character', help='.hic file for condition 2'),
	make_option(c("--cond1"), default=NA, type='character', help='Name of condition 1'),
	make_option(c("--cond2"), default=NA, type='character', help='Name of condition 2'),
	make_option(c("--outdir"), default=NA, type='character', help='Output directory')
)
opt_parser <- OptionParser(option_list = option_list, description = "")
opt <- parse_args(opt_parser)

tad_cond1 <- read.delim(opt$tads1)
tad_cond2 <- read.delim(opt$tads2)
hic_cond1 <- opt$hic_file1
hic_cond2 <- opt$hic_file2
cond1_name <- opt$cond1
cond2_name <- opt$cond2
outDir <- opt$outdir

#tad_cond1 <- read.delim("/mnt/bioadhoc-temp/Groups/vd-ay/dsfigueroa/projects/IKAROS/repository/3D-IKAROS/TAD_calling/WT_TADs_annotated.txt")
#tad_cond2 <- read.delim("/mnt/bioadhoc-temp/Groups/vd-ay/dsfigueroa/projects/IKAROS/repository/3D-IKAROS/TAD_calling/IKDN_TADs_annotated.txt")
#hic_cond1 <- "/mnt/BioAdHoc/Groups/vd-ay/dsfigueroa/projects/IKAROS/datasets/DNAnexus_JUICER_HIC/HIC/NOVASEQ_JUICER_HIC_STDST_mergedRPL/HIC_WT7001_851_Merge_allValidPairsMBOI.hic"
#hic_cond2 <- "/mnt/BioAdHoc/Groups/vd-ay/dsfigueroa/projects/IKAROS/datasets/DNAnexus_JUICER_HIC/HIC/NOVASEQ_JUICER_HIC_STDST_mergedRPL/HIC_IKDN_7002_852_Merge_allValidPairsMBOI.hic"
#cond1_name <- "WT"
#cond2_name <- "IKDN"

resolution <- 10000
matrix_value <- "oe" # oe, observed
diff_thld <- 0

#===============
# Calculate differential score of 
# splitted and merged TADs
#===============

print("Score merged and splitted TADs from condition 1")
scored_tads_cond1 <- Score_Merge_Split_TADs(tad_set1=tad_cond1, tad_set2=tad_cond2, 
											hic_set1=hic_cond1, hic_set2=hic_cond2, 
											cond1=cond1_name, cond2=cond2_name)

splitted_cond1 <- scored_tads_cond1[[1]]
splitted_cond1_filt <- splitted_cond1[splitted_cond1$score_diff > diff_thld,]

merged_cond2 <- scored_tads_cond1[[2]]
merged_cond2 <- merged_cond2[,c(2:4, 1, 5:dim(merged_cond2)[2])]
merged_cond2_filt <- merged_cond2[merged_cond2$score_diff > diff_thld,]

print("Score merged and splitted TADs from condition 2")
scored_tads_cond2 <- Score_Merge_Split_TADs(tad_set1=tad_cond2, tad_set2=tad_cond1, 
											hic_set1=hic_cond2, hic_set2=hic_cond1, 
											cond1=cond2_name, cond2=cond1_name)

splitted_cond2 <- scored_tads_cond2[[1]]
splitted_cond2_filt <- splitted_cond2[splitted_cond2$score_diff > diff_thld,]

merged_cond1 <- scored_tads_cond2[[2]]
merged_cond1 <- merged_cond1[,c(2:4, 1, 5:dim(merged_cond1)[2])]
merged_cond1_filt <- merged_cond1[merged_cond1$score_diff > diff_thld,]

#===============
# Calculate differential score of 
# shifted TADs that share one boundary 
# across conditions
#===============

print("Score shifted TADs")

# Consider a shift of 50kb to calculate score 
window_size <- 50000

shifted_cond1 <- Score_Shift_TADs(tad_set1=tad_cond1, tad_set2=tad_cond2, 
									hic_set1=hic_cond1, hic_set2=hic_cond2, 
									cond1=cond1_name, cond2=cond2_name, slack=window_size)

shifted_cond1_filt <- shifted_cond1[shifted_cond1$score_diff > diff_thld,]

shifted_cond2 <- Score_Shift_TADs(tad_set1=tad_cond2, tad_set2=tad_cond1, 
									hic_set1=hic_cond2, hic_set2=hic_cond1, 
									cond1=cond2_name, cond2=cond1_name, slack=window_size)

shifted_cond2_filt <- shifted_cond2[shifted_cond2$score_diff > diff_thld,]

#===============
# Output scored TADs
#===============

# Remove columns which all values are  NA
merged_cond1 <- merged_cond1[,colSums(is.na(merged_cond1)) < nrow(merged_cond1)]
splitted_cond1 <- splitted_cond1[,colSums(is.na(splitted_cond1)) < nrow(splitted_cond1)]
shifted_cond1 <- shifted_cond1[,colSums(is.na(shifted_cond1)) < nrow(shifted_cond1)]

merged_cond2 <- merged_cond2[,colSums(is.na(merged_cond2)) < nrow(merged_cond2)]
splitted_cond2 <- splitted_cond2[,colSums(is.na(splitted_cond2)) < nrow(splitted_cond2)]
shifted_cond2 <- shifted_cond2[,colSums(is.na(shifted_cond2)) < nrow(shifted_cond2)]

# Save as files
write.table(splitted_cond1, file=paste0(outDir, "/", cond1_name, "_TADs_splitted_in_", cond2_name, "_scored.txt"), 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE) 
write.table(merged_cond1, 
	file=paste0(outDir, "/", cond1_name, "_TADs_merged_in_", cond2_name, "_scored.txt"), 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE) 
write.table(shifted_cond1, file=paste0(outDir, "/", cond1_name, "_TADs_shifted_in_", cond2_name, "_scored.txt"), 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE) 

write.table(splitted_cond2, file=paste0(outDir, "/", cond2_name, "_TADs_splitted_in_", cond1_name, "_scored.txt"), 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE) 
write.table(merged_cond2, file=paste0(outDir, "/", cond2_name, "_TADs_merged_in_", cond1_name, "_scored.txt"), 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE) 
write.table(shifted_cond2, file=paste0(outDir, "/", cond2_name, "_TADs_shifted_in_", cond1_name, "_scored.txt"), 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE) 

#===============
# Perform Aggregate TAD analysis and
# plot Fig3E and SuppFig3E
#===============

print("Perform ATA")

contacts_cond1 <- load_contacts(hic_cond1, 
				resolution = resolution,
				balancing = TRUE,
				sample_name = cond1_name,
				colour = "red")
contacts_cond2 <- load_contacts(hic_cond2, 
				resolution = resolution,
				balancing = TRUE,
				sample_name = cond2_name,
				colour = "red")

# size = 100L for 500kb window in each side at 10kb resolution
window_size <- 100L
contrast_lim <- 2.5

# Plot ATA of unfiltered splitted TADs
ata1_splitted_cond1 <- ATA(contacts_cond1, splitted_cond1, size = window_size)
ata2_splitted_cond1 <- ATA(contacts_cond2, splitted_cond1, size = window_size)

ata1_splitted_cond2 <- ATA(contacts_cond1, splitted_cond2, size = window_size)
ata2_splitted_cond2 <- ATA(contacts_cond2, splitted_cond2, size = window_size)

pdf(paste0(outDir, "/ATA_", cond1_name, "_unfiltered_TADs_splitted_in_", cond2_name, ".pdf"), height=3, width=4)
ContactMapDiff(mat1=ata1_splitted_cond1$signal, mat2=ata2_splitted_cond1$signal, 
		colour_lim=contrast_lim) + labs(title=paste0(cond1_name, " TADs splitted in ", cond2_name))
dev.off()

pdf(paste0(outDir, "/ATA_", cond2_name, "_unfiltered_TADs_splitted_in_", cond1_name, ".pdf"), height=3, width=4)
ContactMapDiff(mat1=ata1_splitted_cond2$signal, mat2=ata2_splitted_cond2$signal, 
		colour_lim=contrast_lim) + labs(title=paste0(cond2_name, " TADs splitted in ", cond1_name))
dev.off()

# Plot ATA of filtered splitted TADs
ata1_splitted_cond1 <- ATA(contacts_cond1, splitted_cond1_filt, size = window_size)
ata2_splitted_cond1 <- ATA(contacts_cond2, splitted_cond1_filt, size = window_size)

pdf(paste0(outDir, "/ATA_", cond1_name, "_filtered_TADs_splitted_in_", cond2_name, ".pdf"), height=3, width=4)
ContactMapDiff(mat1=ata1_splitted_cond1$signal, mat2=ata2_splitted_cond1$signal, 
		colour_lim=contrast_lim) + labs(title=paste0(cond1_name, " TADs splitted in ", cond2_name))
dev.off()

# Plot ATA of merged TADs
ata1_merged_cond1 <- ATA(contacts_cond1, merged_cond1, size = window_size)
ata2_merged_cond1 <- ATA(contacts_cond2, merged_cond1, size = window_size)

ata1_merged_cond2 <- ATA(contacts_cond1, merged_cond2, size = window_size)
ata2_merged_cond2 <- ATA(contacts_cond2, merged_cond2, size = window_size)

pdf(paste0(outDir, "/ATA_", cond1_name, "_TADs_merged_in_", cond2_name, ".pdf"), height=3, width=4)
ContactMapDiff(mat1=ata1_merged_cond1$signal, mat2=ata2_merged_cond1$signal, 
		colour_lim=contrast_lim) + labs(title=paste0(cond1_name, " TADs merged in ", cond2_name))
dev.off()

pdf(paste0(outDir, "/ATA_", cond2_name, "_TADs_merged_in_", cond1_name, ".pdf"), height=3, width=4)
ContactMapDiff(mat1=ata1_merged_cond2$signal, mat2=ata2_merged_cond2$signal, 
		colour_lim=contrast_lim) + labs(title=paste0(cond2_name, " TADs merged in ", cond1_name))
dev.off()

# Plot ATA of contracted TADs in 3' boundary

smaller_tad_cond1 <- shifted_cond1_filt[shifted_cond1_filt$class1==paste0("shifted_smaller_", cond2_name) & 
							shifted_cond1_filt$class2=="shifted_end_boundary",]
smaller_tad_cond2 <- shifted_cond2_filt[shifted_cond2_filt$class1==paste0("shifted_smaller_", cond1_name) & 
							shifted_cond2_filt$class2=="shifted_end_boundary",]

ata1_shifted_cond1 <- ATA(contacts_cond1, smaller_tad_cond1, size = window_size)
ata2_shifted_cond1 <- ATA(contacts_cond2, smaller_tad_cond1, size = window_size)

ata1_shifted_cond2 <- ATA(contacts_cond1, smaller_tad_cond2, size = window_size)
ata2_shifted_cond2 <- ATA(contacts_cond2, smaller_tad_cond2, size = window_size)

pdf(paste0(outDir, "/ATA_", cond1_name, "_TADs_shifted_in_", cond2_name, ".pdf"), height=3, width=4)
ContactMapDiff(mat1=ata1_shifted_cond1$signal, mat2=ata2_shifted_cond1$signal, 
		colour_lim=contrast_lim) + labs(title=paste0(cond1_name, " TADs merged in ", cond2_name))
dev.off()

pdf(paste0(outDir, "/ATA_", cond2_name, "_TADs_shifted_in_", cond1_name, ".pdf"), height=3, width=4)
ContactMapDiff(mat1=ata1_shifted_cond2$signal, mat2=ata2_shifted_cond2$signal, 
		colour_lim=contrast_lim) + labs(title=paste0(cond2_name, " TADs merged in ", cond1_name))
dev.off()

# Plot ATA of filtered out splitted TADs

splitted_cond1_filt_out <- splitted_cond1[splitted_cond1$score_diff <= diff_thld,]

ata1_splitted_cond1 <- ATA(contacts_cond1, splitted_cond1_filt_out, size = window_size)
ata2_splitted_cond1 <- ATA(contacts_cond2, splitted_cond1_filt_out, size = window_size)

pdf(paste0(outDir, "/ATA_", cond1_name, "_filtered_out_TADs_splitted_in_", cond2_name, ".pdf"), height=3, width=4)
ContactMapDiff(mat1=ata1_splitted_cond1$signal, mat2=ata2_splitted_cond1$signal, 
		colour_lim=contrast_lim) + labs(title=paste0(cond1_name, " TADs splitted in ", cond2_name))
dev.off()
