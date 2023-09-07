#!/bin/bash

#===============
# parameters - users need to edit them
#===============

RESOLUTION=50000
H5_MATRIX=sample_${RESOLUTION}.h5
# --minDepth should be at least 3 times as large as the bin size of the Hi-C matrix
MIN_DEPTH=150000 
#--maxDepth should around 6-10 times as large as the bin size of the Hi-C matrix.
MAX_DEPTH=500000 
OUTFILE=sample_${RESOLUTION}_min150k_max500k_step10k_thres0.05_delta0.01_fdr

#===============
# Call topologically associating domains (TADs) 
# using HiCExplorer tool hicFindTADs
#===============

hicFindTADs -m ${H5_MATRIX} --outPrefix ${OUTFILE} --minDepth ${MIN_DEPTH} --maxDepth ${MAX_DEPTH} --step 100000 --thresholdComparisons 0.05 --delta 0.01 --correctForMultipleTesting fdr -p 10

