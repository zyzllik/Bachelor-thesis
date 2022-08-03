#!/bin/sh

# Import inputs from the command line
sample='K562'
path_toTADs="/net/data.isilon/ag-cherrmann/echernova/integrative_nmf/220701_all_windows_filtered.txt.bed"
suffix='_filtered' #whether it is positive or negative, if positive: live blank, if negative add arg '_neg'

path_to_biggies="/net/data.isilon/ag-cherrmann/echernova/bigwig_histone_modifications/$sample"
path_out="/net/data.isilon/ag-cherrmann/echernova/binned_histone_modifications/$sample${suffix}_all-windows"

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

