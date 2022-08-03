#!/bin/sh

n_tads=$1
sample=$2
suffix=$3

tad_size=40000
tad_boundaries="/net/data.isilon/ag-cherrmann/echernova/tad_boundaries_10_40kb/${sample}_TAD_boundaries_central_40kb${suffix}.txt.bed"
output_path="/net/data.isilon/ag-cherrmann/echernova/tad_boundaries_10_40kb/${sample}_negatives_central_40kb${suffix}.txt.bed"

path_dummy="/net/data.isilon/ag-cherrmann/echernova/dummy_bed.txt.bed"

counter=0

rm $output_path

for i in $(seq 1 $n_tads)
do
((counter++))
echo -e "chr1\t1\t40001" >> $path_dummy
done

echo $counter

source ~/miniconda3/bin/activate
conda activate testenv

bedtools shuffle -i $path_dummy -g hg38_imputed.genome -excl $tad_boundaries >> $output_path

rm $path_dummy
