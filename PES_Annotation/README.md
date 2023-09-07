Promoters (P), enhancer (E) and structural (S) loop anchor annotation
=============================

Loop anchor annotation into promoters (P), enhancer (E) and structural (S) given 
a set of interactions, CTCF, H3K27ac and H3K4me1 ChIP-seq.

Usage: PES_Annotation.R [options]

Options:

	--loops
		Input file with loops coordinates in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format

	--genes
		File with genes TSS coordinates

	--ctcf
		File with CTCF ChIP-seq peak coordinates

	--h3k27ac
		File with H3K27ac ChIP-seq peak coordinates

	--h3k4me1
		File with H3K4me1 ChIP-seq peak coordinates
	
	--slack
		bp of shift allowed for overlap

	--outfile
		Output filename

