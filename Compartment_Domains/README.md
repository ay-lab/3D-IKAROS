# Step-by-step execution

## 1. Compartment domains

Define compartment domains diven a set of subcompartments.

Subcompartments calling was done with a modified functions from CALDER v1.0 to take mm10 genome.Cuurently CALDER 2 support multiple genomeversions. 

Usage: Define_CompartmentDomains.R [options]

Options:

	--subcomp
		Subcompartment outfile from CALDER (*_sub_compartments.bed)
	
	--sample_id
		Sample ID

	--outdir
		Output directory


## 2.1 Intra- or inter-compartmental interactions

Classify interactions into intra- or inter-compartmental interactions

Usage: Overlap_CompartmentDomains_Loops.R [options]

Options:

	--domain_cond1
		Compartment domains in condition 1. Outfile from ./Define_CompartmentDomains.R (*_CompartmentDomains.bed)
	
	--domain_cond2
		Compartment domains in condition 2. Outfile from ./Define_CompartmentDomains.R (*_CompartmentDomains.bed)

	--cond1
		Name of condition 1

	--cond2
		Name of condition 2
	
	--up_loops
		File with upregulated loops coordinates in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format

	--down_loops
		File with downregulated loops coordinates in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format

	--resolution
		Bin size

	--outdir
		Output directory


## 2.2 Compartmental change at loop anchors

Compute the change of PC values of compartments overlapping loop anchors from differential interactions

Usage: Compartmental_changes_Diff_Intercations.R [options]

Options:

	--comp
		Compartment from dcHiC output (intra_compartment.bedGraph)
	
	--cond1
		Name of condition 1

	--cond2
		Name of condition 2
	
	--up_loops
		File with upregulated loops coordinates in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format

	--down_loops
		File with downregulated loops coordinates in "chr1 \t start1 \t end1 \t chr2 \t start2 \t end2" format

	--outdir
		Output directory


