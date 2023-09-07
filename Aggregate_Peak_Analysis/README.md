Aggregate Peak Analysis (APA)
=============================

Compute aggregate enrichment of a set of interactions in a contact matrix

Defined in the paper:
Rao et. al., A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping, Cell 2014.

Usage: APA_GENOVA.R [options]

Options:

	--loops
		Input file with loops coordinates in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format

	--hic_file
		Contact map in .hic format

	--chr_name
		If 1, the chromosome names in the contact map are in "chr#" format. If 0, they are in "#" format. Default 1

	--min_lim
		Lower limit of color scale

	--max_lim
		Upper limit of color scale
	
	--outfile
		Output filename

