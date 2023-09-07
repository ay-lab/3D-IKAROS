Overlap of genes and loops
=============================

Usage: Overlap_genes_loops.R [options]

Options:

	--genes
		File with genes TSS coordinates

	--dgea_file
		Differential Gene Expression Analysis output file with "gene_name \t adj_pvalue \t log2_fold_change" format

	--up_loops
		File with upregulated loops coordinates in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format

	--down_loops
		File with downregulated loops coordinates in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format

	--nd_loops
		File with non-differential loops coordinates in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format

	--outdir
		Output directory

