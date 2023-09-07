# Step-by-step execution

## 1. Convert HiCPro output into H5 format

./convert_HiCPro_format_into_h5.sh is for converting HiCPro output into H5 format using hicConvertFormat tool from HiCExplorer

In the parameter section, specify:

MATRIX: ICE normalized contact map
BEDFILE: BED file of genomic bins
RESOLUTION: Bin size
OUTFILE: Filename of h5 matrix

## 2. TAD calling

./call_TADs.sh is for calling topologically associating domains (TADs) using hicFindTADs tool from HiCExplorer

In the parameter section, specify:

H5_MATRIX: Contact map in h5 format
RESOLUTION: Bin size
MIN_DEPTH: Minimum window length. This number should be at least 3 times as large as the bin size of the Hi-C matrix
MAX_DEPTH: Maximum window length. This number should around 6-10 times as large as the bin size of the Hi-C matrix. 
OUTFILE: Filename ID

## 3. TAD classification

Classify TADs into common, splitted, merged and shifted

Usage: TAD_classification.R [options]

Options:

	--tads1
		File with TAD coordinates from condition 1 with "chr \t start \t end" format

	--tads2
		File with TAD coordinates from condition 2 with "chr \t start \t end" format

	--cond1
		Name of condition 1

	--cond2
		Name of condition 2

	--outdir
		Output directory

## 4. TAD_filtering.R

Filter TADs based on differences in contacts density between the two conditions and perform aggregate TAD analysis

Usage: TAD_filtering.R [options]

Options:

	--tads1
		File with TAD coordinates from condition 1 with "chr \t start \t end" format

	--tads2
		File with TAD coordinates from condition 2 with "chr \t start \t end" format

    --hic_file1
		Contact map from condition 1 in .hic format

	--hic_file2
		Contact map from condition 2 in .hic format

	--cond1
		Name of condition 1

	--cond2
		Name of condition 2

	--outdir
		Output directory

