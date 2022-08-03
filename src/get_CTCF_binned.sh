#!/bin/sh

## This script will:
## - make individual summarised versions of the windows defined as TAD boundaries using epigenomic features (and RNA-seq signal)
## - we can then bind the windows together and add in the results from the CTCF footprinting

# We will use the same bigwigs processed uniformly that we used for the NMF analysis
path_to_biggies=/net/data.isilon/ag-cherrmann/echernova/bigwig_histone_modifications/GM12878
path_toTADs=/net/data.isilon/ag-cherrmann/echernova/tad_boundaries_10_40kb/GM12878_TAD_boundaries_10_40kb.txt.bed
path_out=/net/data.isilon/ag-cherrmann/echernova/binned_histone_modifications/GM12878

bw=$(find $path_to_biggies -name "*CTCF*")
echo $bw
base="${bw%*.bigWig}"
label=${base##*/}
tadfile_id="${path_toTADs%*.bed}"
tadf=${tadfile_id##*/}

source ~/miniconda3/bin/activate
conda activate testenv
multiBigwigSummary BED-file -b $bw -out $path_out/${label}_${tadf}.npz --labels ${label%*_RPKM} --outRawCounts $path_out/${label}_${tadf}.tab --BED $path_toTADs

