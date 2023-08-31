#!/bin/bash

#===============
# Convert HiCPro output into H5 format 
# using HiCExplorer tool hicConvertFormat
#===============

MATRIX=sample_iced.matrix
BEDFILE=sample_abs.bed
RESOLUTION=50000
OUTFILE=sample_${RESOLUTION}.h5

hicConvertFormat -m ${MATRIX} --inputFormat hicpro --bedFileHicpro ${BEDFILE} --outputFormat h5 --outFileName ${OUTFILE}
