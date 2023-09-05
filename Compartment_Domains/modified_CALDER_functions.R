#!/usr/bin/env Rscript

#===========================================================
# The function 'get_gene_info' from the R package CALDER v1.0 
# was modified to take mm10 genome build as input

# Author: Daniela Salgado Figueroa
#===========================================================

#================================================
CALDER_CD_hierarchy_v2 = function(contact_mat_file, 
							   chr, 
							   bin_size, 
							   out_dir, 
							   save_intermediate_data=FALSE,
							   genome = 'hg19')
{
    time0 = Sys.time()
    log_file = paste0(out_dir, '/chr', chr, '_log.txt')

    cat('\n')

    cat('>>>> Begin process contact matrix and compute correlation score at:', as.character(Sys.time()), '\n', file=log_file, append=FALSE)
    cat('>>>> Begin process contact matrix and compute correlation score at:', as.character(Sys.time()), '\n')
    processed_data = contact_mat_processing(contact_mat_file, bin_size=bin_size)
   
    A_oe = processed_data$A_oe
    ccA_oe_compressed_log_atanh = processed_data$atanh_score

    cat('\r', '>>>> Finish process contact matrix and compute correlation score at:', as.character(Sys.time()))
    cat('>>>> Finish process contact matrix and compute correlation score at:', as.character(Sys.time()), '\n', file=log_file, append=TRUE)

    p_thresh = ifelse(bin_size < 40000, 0.05, 1)
    window.sizes = 3
    compartments = vector("list", 2)
    chr_name = paste0("chr", chr)

    cat('>>>> Begin compute compartment domains and their hierachy at:', as.character(Sys.time()), '\n', file=log_file, append=TRUE)
    cat('\r', '>>>> Begin compute compartment domains and their hierachy at:', as.character(Sys.time()))

    compartments[[2]] = generate_compartments_bed(chr = chr, bin_size = bin_size, window.sizes = window.sizes, ccA_oe_compressed_log_atanh, p_thresh, out_file_name = NULL, stat_window_size = NULL)
    topDom_output = compartments[[2]]
    bin_names = rownames(A_oe)
    A_oe = as.matrix(A_oe)
    initial_clusters = apply(topDom_output$domain[, c("from.id", "to.id")], 1, function(v) v[1]:v[2])

    if (sum(sapply(initial_clusters, length)) != max(unlist(initial_clusters))) {
        stop(CELL_LINE, " initial_clusters error in topDom")
    }

    n_clusters = length(initial_clusters)
		A_oe_cluster_mean = HighResolution2Low_k_rectangle(A_oe, initial_clusters, initial_clusters, sum_or_mean = "mean")

	trend_mean_list = lapply( 1:4, function(v) 1*(A_oe_cluster_mean[, -(1:v)] > A_oe_cluster_mean[, - n_clusters - 1 + (v:1)]) )
	trend_mean = do.call(cbind, trend_mean_list)
	c_trend_mean = cor(t(trend_mean))
	atanh_c_trend_mean= atanh(c_trend_mean / (1+1E-7))


	# if(to_scale)
	{
		trend_mean = scale(trend_mean)
		c_trend_mean = scale(c_trend_mean)
		atanh_c_trend_mean= scale(atanh_c_trend_mean)
	}


	PC_12_atanh = get_PCs(atanh_c_trend_mean, which=1:10)
	PC_12_atanh[, 2:10] = PC_12_atanh[, 2:10]/5 ## xx-xx-xxxx: compress PC2
	rownames(PC_12_atanh) = 1:nrow(PC_12_atanh)

	############################################################
	PC_direction = 1

	gene_info <- get_gene_info_v2(genome)

	## switch PC direction based on gene density
	{
		initial_clusters_ori_bins = lapply(initial_clusters, function(v) as.numeric(bin_names[v]))
		chr_bin_pc = data.table::data.table(chr=chr_name, bin=unlist(initial_clusters_ori_bins), PC1_val=rep(PC_12_atanh[,1], sapply(initial_clusters_ori_bins, length)))
		chr_bin_pc$start = (chr_bin_pc$bin - 1)*bin_size + 1
		chr_bin_pc$end = chr_bin_pc$bin*bin_size
		chr_bin_pc_range = makeGRangesFromDataFrame(chr_bin_pc, keep.extra.columns=TRUE)
		gene_info_chr = subset(gene_info, seqnames==chr_name)

		refGR = chr_bin_pc_range
		testGR = gene_info_chr
		hits <- findOverlaps(refGR, testGR)
	    overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
	    overlaps_bins = unique(data.table::data.table(overlap_ratio=width(overlaps)/bin_size, bin=overlaps$bin))
	    bin_pc_gene_coverage = merge(chr_bin_pc, overlaps_bins, all.x=TRUE)
	    bin_pc_gene_coverage$overlap_ratio[is.na(bin_pc_gene_coverage$overlap_ratio)] = 0
		
		gene_density_cor = cor(method='spearman', subset(bin_pc_gene_coverage, (PC1_val < quantile(PC1_val, 0.25)) | (PC1_val > quantile(PC1_val, 0.75)) , c('PC1_val', 'overlap_ratio')))[1,2]
	    
	    if(abs(gene_density_cor) < 0.2) warning('correlation between gene density and PC1 is too weak')
	    PC_direction = PC_direction*sign(gene_density_cor)

	    PC_12_atanh = PC_12_atanh*PC_direction
	}


	project_info = project_to_major_axis(PC_12_atanh)
	x_pro = project_info$x_pro
	
	############################################################
	hc_disect_kmeans_PC12 = bisecting_kmeans(PC_12_atanh[, 1:10, drop=FALSE])

	hc_hybrid_PC12 = hc_disect_kmeans_PC12

	{
		reordered_names = reorder_dendro(hc_hybrid_PC12, x_pro, aggregateFun=mean)
		hc_hybrid_PC12_reordered = dendextend::rotate(hc_hybrid_PC12, reordered_names)
		hc_hybrid_x_pro = hc_disect_kmeans_PC12
		reordered_names_x_pro = get_best_reorder(hc_hybrid_x_pro, x_pro)
		CALDER_hc = dendextend::rotate(hc_hybrid_x_pro, reordered_names_x_pro)	
	}

	############################################################
	parameters = list(bin_size = bin_size, p_thresh = p_thresh)
	res = list( CALDER_hc=CALDER_hc, initial_clusters=initial_clusters, bin_names=bin_names, x_pro=x_pro, parameters=parameters)
	print(res)
	intermediate_data_file = paste0(out_dir, '/chr', chr, '_intermediate_data.Rds')
	
	hc = res$CALDER_hc
	hc_k_labels_full = try(get_cluser_levels(hc, k_clusters=Inf, balanced_4_clusters=FALSE)$cluster_labels)
	bin_comp = data.table::data.table(chr=chr, bin_index=res$bin_names, comp=rep(hc_k_labels_full, sapply(res$initial_clusters, length)))

	rownames(bin_comp) = NULL
	res$comp = bin_comp
	res$CDs = lapply(res$initial_clusters, function(v) res$bin_names[v])
	res$mat = A_oe
	res$chr = chr
	generate_hierachy_bed(chr=chr, res=res, out_dir=out_dir, bin_size=bin_size)


    cat('>>>> Finish compute compartment domains and their hierachy at: ', as.character(Sys.time()), '\n', file=log_file, append=TRUE)
    cat('\r', '>>>> Finish compute compartment domains and their hierachy at: ', as.character(Sys.time()))

   	if(abs(gene_density_cor) < 0.2) cat('WARNING: correlation between gene density and PC1 on this chr is: ', substr(gene_density_cor, 1, 4), ', which is a bit low', '\n', file=log_file, append=TRUE)

    time1 = Sys.time()
    # delta_time  = gsub('Time difference of', 'Total time used for computing compartment domains and their hierachy:', print(time1 - time0))

    delta_time <- time1 - time0
	timediff <- format(round(delta_time, 2), nsmall = 2)

    cat('\n\n', 'Total time used for computing compartment domains and their hierachy:', timediff, '\n', file=log_file, append=TRUE)
   	# if(abs(gene_density_cor) > 0.2) cat('The gene density and PC1 correlation on this chr is: ', substr(gene_density_cor, 1, 4), '\n', file=log_file, append=TRUE)

	############################################################
	intermediate_data = res
	if(save_intermediate_data==TRUE) saveRDS(intermediate_data, file=intermediate_data_file)
	# cat(intermediate_data_file)
	return(intermediate_data)
}

get_gene_info_v2 <- function(genome)
{

	if(genome == 'hg19') {
		gene_info_file = system.file("extdata", "TxDb.Hsapiens.UCSC.hg19.knownGene.rds", package = 'CALDER')
	} else if(genome == 'mm9') {
		gene_info_file = system.file("extdata", "TxDb.Mmusculus.UCSC.mm9.knownGene.rds", package = 'CALDER')
	}else if(genome == 'mm10') {
		gene_info_file = system.file("extdata", "TxDb.Mmusculus.UCSC.mm10.knownGene.rds", package = 'CALDER')
	} else {
		print("Hello")
		stop(paste0("Unknown genome (", genome, ")"))

	}

	
	gene_info = readRDS(gene_info_file)
}

CALDER_sub_domains_v2 = function(intermediate_data_file=NULL, intermediate_data=NULL, chr, out_dir, bin_size)
{	
    time0 = Sys.time()
    log_file = paste0(out_dir, '/chr', chr, '_sub_domains_log.txt')

   	cat('\r', '>>>> Begin compute sub-domains at:', as.character(Sys.time()))
   	cat('>>>> Begin compute sub-domains at:', as.character(Sys.time()), '\n', file=log_file, append=FALSE)

	if(is.null(intermediate_data)) intermediate_data = readRDS(intermediate_data_file)
	{
	    if(intermediate_data$chr!=chr) stop('intermediate_data$chr!=chr; check your input parameters\n') 
	    if( !setequal(rownames(intermediate_data$mat), intermediate_data$bin_names) ) stop('!setequal(rownames(intermediate_data$mat), intermediate_data$bin_names) \n')     
	    compartment_segs = generate_compartment_segs( intermediate_data$initial_clusters )

			cat('\r', '>>>> Begin compute sub-domains within each compartment domain at:', as.character(Sys.time()))   			
		cat('>>>> Begin compute sub-domains within each compartment domain at:', as.character(Sys.time()), '\n', file=log_file, append=TRUE)
		sub_domains_raw = HRG_zigzag_compartment_domain_main_fun(intermediate_data$mat, './', compartment_segs, min_n_bins=2)   
	    no_output = post_process_sub_domains(chr, sub_domains_raw, ncores=1, out_dir=out_dir, bin_size=bin_size)
	    cat('>>>> Finish compute sub-domains within each compartment domain at:', as.character(Sys.time()), '\n', file=log_file, append=TRUE)
	    cat('\r', '>>>> Finish compute sub-domains within each compartment domain at:', as.character(Sys.time()))

	   	time1 = Sys.time()
        # delta_time  = gsub('Time difference of', 'Total time used for computing compartment domains and their hierachy:', print(time1 - time0))
        delta_time <- time1 - time0
		timediff <- format(round(delta_time, 2), nsmall = 2)

        cat('\n\n', 'Total time used for computing sub-domains:', timediff, '\n', file=log_file, append=TRUE)
	}
	# return(NULL)
}

CALDER_main_v2 = function(contact_mat_file, 
					   chr, 
					   bin_size, 
					   out_dir, 
					   sub_domains = TRUE, 
					   save_intermediate_data = TRUE,
					   genome='hg19') {
	required_packages = c('doParallel', 'GenomicRanges', 'R.utils', 'factoextra', 'maptools')
	sapply(required_packages, require, character.only = TRUE, quietly = TRUE)

	dir.create(out_dir, showWarnings = FALSE)
	intermediate_data = CALDER_CD_hierarchy_v2(contact_mat_file, 
											chr, 
											bin_size, 
											out_dir, 
											save_intermediate_data,
											genome)

	if(sub_domains==TRUE) {
		CALDER_sub_domains_v2(intermediate_data=intermediate_data, 
						   chr=chr, 
						   out_dir=out_dir, 
						   bin_size=bin_size)
	}
}

#================================================

# Add genome to package directory

library(TxDb.Mmusculus.UCSC.mm10.knownGene)

package_dir <- "/mnt/BioHome/dsfigueroa/micromamba/envs/R-4.1.3/lib/R/library/CALDER"

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)

saveRDS(genes, file=paste0(package_dir, "/extdata/TxDb.Mmusculus.UCSC.mm10.knownGene.rds"))
