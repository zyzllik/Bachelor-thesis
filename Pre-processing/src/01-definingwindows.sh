#!/bin/sh

## This script will:
## - make individual summarised versions of the windows defined as TAD boundaries using epigenomic features (and RNA-seq signal)
## - we can then bind the windows together and add in the results from the CTCF footprinting

# We will use the same bigwigs processed uniformly that we used for the NMF analysis
# path_to_biggies=/net/data.isilon/ag-cherrmann/echernova/bigwig_histone_modifications/IMR90
# path_toTADs=/net/data.isilon/ag-cherrmann/echernova/tad_boundaries_10_40kb/IMR90_TAD_boundaries_10_40kb.txt.bed
# path_out=/net/data.isilon/ag-cherrmann/echernova/binned_histone_modifications/IMR90

# Import inputs from the command line
sample=$1
path_toTADs=$2
suffix=$3 #whether it is positive or negative, if positive: live blank, if negative add arg '_neg'

path_to_biggies="/net/data.isilon/ag-cherrmann/echernova/bigwig_histone_imputed/$sample"
path_out="/net/data.isilon/ag-cherrmann/echernova/binned_histone_modifications/$sample$suffix"

mkdir $path_out

# the loop for submission is ready:
for bw in $path_to_biggies/*.bigWig
do 
base="${bw%*.bigWig}"
label=${base##*/}
tadfile_id="${path_toTADs%*.bed}"
tadf=${tadfile_id##*/}

sbatch --wrap="source ~/miniconda3/bin/activate; \
conda activate testenv; \
multiBigwigSummary BED-file -b $bw -out $path_out/${label}_${tadf}.npz --labels ${label%*_RPKM} --outRawCounts $path_out/${label}_${tadf}.tab --BED $path_toTADs "
done

