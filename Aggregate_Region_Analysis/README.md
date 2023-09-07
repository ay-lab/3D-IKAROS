Aggregate Region Analysis (ARA)
=============================

Aggregate enrichment of a set of regions in a contact matrix

Usage: ARA_GENOVA.R [options]

Options:

	--regions
		Input file with regions coordinates in "chr \t start \t end" format

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
	