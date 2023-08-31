## Step-by-step execution

1. convert_HiCPro_format_into_h5.sh - Convert HiCPro output into H5 format using hicConvertFormat tool from HiCExplorer
2. call_TADs.sh - Call topologically associating domains (TADs) using hicFindTADs tool from HiCExplorer
3. TAD_classification.R - Classify TADs into common, splitted, merged and shifted
4. TAD_filtering.R - Filter TADs based on differences in contacts density between the two conditions and perform aggregate TAD analysis

