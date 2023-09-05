#!/usr/bin/env Rscript

#===========================================================
# R script for performing a time course analysis
# given a set of interactions. The following code was copied 
# from Tcseq package and modified by Abhijit Chakraborty & 
# Daniela Salgado Figueroa
# Wu M, Gu L (2021). TCseq: Time course sequencing data analysis. R package version 1.18.0.
#===========================================================

#================================================
clust <- methods::setClass("clust", slots = c(method = "character",
				dist = "character",
				data = "matrix",
				centers = "matrix",
				cluster = "integer",
				membership = "matrix"))

pcaclust <- function(x, algo, k, dist = "euclidean", centers = NULL,
	standardize = TRUE, ...) {
	if (is.matrix(x)) {
	data.tmp <- x
	}else{
	data.tmp <- x@tcTable
	}
	if (standardize) {
	for (i in seq_len(nrow(data.tmp))) {
	data.tmp[i, ] <- (data.tmp[i, ] - mean(data.tmp[i, ], na.rm = TRUE))/sd(data.tmp[i, ], na.rm = TRUE)
	}
	data.tmp <- data.tmp[complete.cases(data.tmp), ]
	}
	object <- methods::new("clust")
	object@method <- algo
	object@dist <- dist
	object@data <- data.tmp
	
	res <- .pcaclust(data = data.tmp, algo = algo, k = k, dist = dist,
	centers = centers, ...)
	
	if (algo == "cm") {
	object@cluster <- res$cluster
	object@membership <- res$membership
	object@centers <- res$centers
	} else {
	object@cluster <- res$cluster
	object@centers <- res$centers
	}
	if (is.matrix(x)) {
	object
	} else {
	x@clusterRes <- object
	x
	}
}

# perform time course clustering
.pcaclust <- function(data, algo, k, centers = NULL,
	dist = "euclidean", ...) {
	if (!algo %in% c("pam", "km", "hc", "cm")) {
	stop("clustering method should be one of 'pam','km','hc','cm'")
	}
	if (!dist %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
	stop("Distance metric should be 'correlation', or one of the distance measures in dist function")
	}
	if (algo == "km") {
	if(dist != "euclidean"){
	stop("kmeans only support euclidean metric; for other distance metrices, please see the help page")
	}
	}
	if (algo == "cm" ) {
	if(!dist %in% c("euclidean", "manhattan")){
	stop("cmeans only support euclidean or mahattan distance metrics")
	}
	}
	
	d <- NULL
	if (algo %in% c("pam", "hc")) {
	if (dist == "correlation") {
	d <- as.dist(1 - cor(t(data)))
	}
	if (dist != "correlation") {
	d <- dist(data, method = dist)
	}
	}
	clustres <- list()
	if (algo != "hc") {
	if (!is.null(centers)) {
	if (nrow(centers) != k) {
	stop("Number of rows of centers must be equal to k")
	}
	}
	}
	clustres <- switch(algo, km = {
	if (!is.null(centers)) {
	res <- kmeans(data, centers = centers, ...)
	} else {
	res <- kmeans(data, centers = k, ...)
	}
	clustres$cluster <- res$cluster
	clustres$centers <- res$centers
	clustres
	}, pam = {
	if (!is.null(centers)) {
	ind <- data[, 1] %in% centers[, 1]
	ind <- which(ind)
	if (length(ind) != k) {
	stop("For 'pam', centers must be chosen from the data")
	} else {
	res <- pam(d, k = k, medoids = ind, ...)
	}
	}
	res <- pam(d, k = k, ...)
	clustres$cluster <- res$clustering
	clustres$centers <- data[res$medoids, ]
	clustres
	}, hc = {
	tree <- hclust(d, ...)
	res <- cutree(tree, k = k)
	clustres$cluster <- res
	clustres$centers <- matrix(0, 0, 0)
	clustres
	}, cm = {
	if (!is.null(centers)) {
	res <- e1071::cmeans(data, centers = centers, ...)
	} else {
	res <- e1071::cmeans(data, centers = k, ...)
	}
	clustres$cluster <- res$cluster
	clustres$centers <- res$centers
	clustres$membership <- res$membership
	clustres
	})
	clustres
}

pcaclustplot <- function(object = NULL, categories = "timepoint",
	value = "expression", cols = NULL,
	cl.color = "gray50",
	membership.color = rainbow(30, s = 3/4, v = 1, start = 1/6),
	title.size = 18, axis.line.size = 0.6,
	axis.title.size = 18,
	axis.text.size = 16, legend.title.size = 14,
	legend.text.size = 14, outprefix= "pca_clusters") {

	if (class(object) != "clust" && class(object) != "TCA") {
	stop("object should be a 'timeclust' object or a 'TCA' object")
	}
	if (class(object) == "clust") {
	data <- object@data
	cluster <- object@cluster
	membership <- object@membership
	}
	if (class(object) == "TCA") {
	data <- object@clusterRes@data
	cluster <- object@clusterRes@cluster
	membership <- object@clusterRes@membership
	}
	ncl <- max(cluster)
	membercolor <- vector(length = length(cluster))
	membervalue <- list()
	counter <- 0
	if (!sum(dim(membership) == 0) == 2) {
	color <- membership.color
	colorseq <- seq(0, 1, length = length(color))
	for (i in seq_len(ncl)) {
	mtmp <- membership[cluster == i, i]
	membervalue[[i]] <- mtmp
	for (j in seq_len(length(mtmp))) {
	counter <- counter + 1
	ind <- which(abs(colorseq - mtmp[j]) == min(abs(colorseq - mtmp[j])))
	membercolor[counter] <- color[ind]
	}
	}
	membervalue <- unlist(membervalue)
	names(membercolor) <- membervalue
	}

	plotlist <- list()
	for (i in seq_len(ncl)) {
	title <- paste0("Cluster ", i)
	dtmp <- data[cluster == i, ]
	a <- which(cluster == i)
	if (length(a) == 1) {
	dtmp <- data.frame(time = 1:length(dtmp), value = dtmp)
	if (!sum(dim(membership) == 0) == 2) {
	m <- membership[cluster == i, i]
	colorname = toString(m)
	plotlist[[i]] <- ggplot2::ggplot(dtmp, ggplot2::aes(x = time, y = value)) +
	ggplot2::geom_line(colour = membercolor[colorname]) + ggplot2::theme_bw() +
	ggplot2::ggtitle(title) +
	ggplot2::scale_x_continuous(breaks = dtmp$time,
	labels = row.names(dtmp)) +
	ggplot2::labs(x = categories, y = value) +
	ggplot2::theme(plot.title = ggplot2::element_text(size = title.size),
	axis.line.x = ggplot2::element_line(color = "black",
	size = axis.line.size),
	axis.line.y = ggplot2::element_line(color = "black",
	size = axis.line.size),
	axis.title = ggplot2::element_text(size = axis.title.size),
	axis.text = ggplot2::element_text(size = axis.text.size),
	legend.position = "none", panel.border = ggplot2::element_blank(),
	panel.grid.major = ggplot2::element_blank(),
	panel.grid.minor = ggplot2::element_blank())
	} else {
	plotlist[[i]] <- ggplot2::ggplot(dtmp, ggplot2::aes(x = time, y = value)) +
	ggplot2::geom_line(colour = cl.color) + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
	ggplot2::scale_x_continuous(breaks = dtmp$time,
	labels = row.names(dtmp)) +
	ggplot2::labs(x = categories, y = value) +
	ggplot2::theme(plot.title = ggplot2::element_text(size = title.size),
	axis.line.x = ggplot2::element_line(color = "black",
	size = axis.line.size),
	axis.line.y = ggplot2::element_line(color = "black",
	size = axis.line.size),
	axis.title = ggplot2::element_text(size = axis.title.size),
	axis.text = ggplot2::element_text(size = axis.text.size),
	legend.position = "none", panel.border = ggplot2::element_blank(),
	panel.grid.major = ggplot2::element_blank(),
	panel.grid.minor = ggplot2::element_blank())
	}
	} else {
	dtmp_m <- reshape2::melt(dtmp)
	colnames(dtmp_m) <- c("group", "time", "value")
	if (sum(dim(membership) == 0) == 2) {
	plotlist[[i]] <- ggplot2::ggplot(dtmp_m, ggplot2::aes(x = time, y = value)) +
	ggplot2::geom_line(ggplot2::aes(group = group), colour = cl.color) +
	ggplot2::theme_bw() + ggplot2::ggtitle(title) +
	ggplot2::labs(x = categories, y = value) +
	ggplot2::theme(plot.title = ggplot2::element_text(size = title.size),
	axis.line.x = ggplot2::element_line(color = "black",
	size = axis.line.size),
	axis.line.y = ggplot2::element_line(color = "black",
	size = axis.line.size),
	axis.title = ggplot2::element_text(size = axis.title.size),
	axis.text = ggplot2::element_text(size = axis.text.size),
	legend.position = "none", panel.border = ggplot2::element_blank(),
	panel.grid.major = ggplot2::element_blank(),
	panel.grid.minor = ggplot2::element_blank())
	}
	if (!sum(dim(membership) == 0) == 2) {
	mem <- membership[cluster == i, i]
	mem1 <- data.frame(group = names(mem), member = mem)
	dtmp_m1 <- merge(dtmp_m, mem1, by = "group")
	colnames(dtmp_m1) <- c("group", "time", "value", "membership")
	dtmp_m1 <- dtmp_m1[order(dtmp_m1[, 4]), ]
	new.factor <- unique(as.vector(dtmp_m1$group))
	dtmp_m1$group <- factor(dtmp_m1$group, levels = new.factor)

	plotlist[[i]] <- ggplot2::ggplot(dtmp_m1, ggplot2::aes(x = time, y = value,
	colour = membership)) +
	ggplot2::geom_line(ggplot2::aes(group = group)) +
	ggplot2::scale_colour_gradientn(colours = membership.color) +
	ggplot2::guides(colour = ggplot2::guide_colourbar()) + ggplot2::theme_bw() +
	ggplot2::ggtitle(title) + ggplot2::labs(x = categories, y = value) +
	ggplot2::theme(plot.title = ggplot2::element_text(size = title.size),
	axis.line.x = ggplot2::element_line(color = "black",
	size = axis.line.size),
	axis.line.y = ggplot2::element_line(color = "black",
	size = axis.line.size),
	axis.title = ggplot2::element_text(size = axis.title.size),
	axis.text = ggplot2::element_text(size = axis.text.size),
	legend.title = ggplot2::element_text(size = legend.title.size),
	legend.text = ggplot2::element_text(size = legend.title.size),
	panel.border = ggplot2::element_blank(),
	panel.grid.major = ggplot2::element_blank(),
	panel.grid.minor = ggplot2::element_blank())


	}
	}
	}
	plots <- c(list, plotlist)
	for(i in 1:length(plots)) {
	print(plots[[i]])
	}
	dev.off() 
}

#================================================

#===============
# Read loop files
#===============

k27ac_up <- read.delim("Data/Upregulated_EdgeR_glmQLFTest_H3K27ac_loops.txt")
k27ac_down <- read.delim("Data/Downregulated_EdgeR_glmQLFTest_H3K27ac_loops.txt")

ctcf_up <- read.delim("Data/Upregulated_EdgeR_glmQLFTest_CTCF_loops.txt")
ctcf_down <- read.delim("Data/Downregulated_EdgeR_glmQLFTest_CTCF_loops.txt")

outDir <- "./"

#===============
# Cluster loops based on 
# observed/expected (OE) values across
# days of IKAROS induction
#===============

interactions <- list(Up_H3K27AC=k27ac_up, Down_H3K27AC=k27ac_down, 
				Up_CTCF=ctcf_up, Down_CTCF=ctcf_down)
kcenter_list <- c(5,5,6,6)
minprob <- 0

set.seed(5)
for(i in 1:length(interactions)){
	
	loops_df <- interactions[[i]]
	pcafile <- loops_df
	pcafile[is.na(pcafile)] <- 0
	
	kcenter <- kcenter_list[i]

	loop_id <- names(interactions)[i]
	print(loop_id)

	# Extract O/E values across days of IKAROS induction
	samples <- grep("D0|D3|D12", colnames(loops_df), value=TRUE)
	samples_id <- c("D0", "D3", "D12")

	pcafile <- pcafile[,c("chr1", "start1", "start2", "loop_id", samples)]
	pcafile <- pcafile[order(pcafile$chr1, pcafile$start1),]
	rownames(pcafile) <- pcafile$loop_id
	pcafile <- as.matrix(pcafile[,-c(1:4)])
	colnames(pcafile) <- samples_id
	head(pcafile)

	# Perform soft k-means clustering
	pcacl <- pcaclust(pcafile, algo="cm", k=kcenter, standardize = TRUE)

	# Plot loop clusters
	plotname <- paste0(outDir, "/TCseq_", loop_id, "_clusters.pdf") 
	pcaclustplot(pcacl, value="Z-score( O/E )", cols = 1)
	system(paste0("mv Rplots.pdf ", plotname), wait = TRUE)

	# Filter loops by membership (?)
	cl_df <- as.data.frame(pcacl@membership)
	cl_df <- cl_df[apply(cl_df, 1, max) > minprob,]
	cat("Keeing members that are above minimum probability of ", minprob, "\n")
	cl_df[,"cluster"] <- paste0("Cluster_",apply(cl_df, 1, which.max))
	cl_df[,"loop_id"] <- rownames(cl_df)
	
	table(cl_df$cluster)

	# Save table of loop with cluster assignation
	cols <- c("loop_id", "cluster", 1:kcenter)
	filename <- paste0(outDir, "/TCseq_", loop_id, "_clusters.txt") 
	
	write.table(cl_df[,cols], file = filename, quote = FALSE, sep = "\t", 
	row.names = FALSE, col.names = TRUE) 

}


