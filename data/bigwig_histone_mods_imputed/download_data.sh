#!/bin/bash
cell_ids=(E114 E116 E118 E028 E122 E017 E123 E127)
cell_names=(A549 GM12878 HepG2 HMEC HUVEC IMR90 K562 NHEK)
histone_mods=(
    "H3K4me1"
    "H3K4me3"
    "H3K36me3"
    "H3K27me3"
    "H3K9me3"
    "H3K27ac"
    "H3K9ac"
    "H3K4me2"
    "H3K79me2"
    "H4K20me1"
    'H2A.Z'
)

n_samples=${#cell_ids[@]}

for ((i=0; i<$n_samples; i++))
do
    cell_n=${cell_names[$i]}
    c_id=${cell_ids[$i]}

    # make a new directory for each cell and enter it
    mkdir $cell_n
    cd $cell_n

    # download data
    for hist_m in ${histone_mods[@]}
    do
        url="https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/$hist_m/$c_id-$hist_m.imputed.pval.signal.bigwig"
        wget --no-check-certificate $url
        mv "$c_id-$hist_m.imputed.pval.signal.bigwig" "$cell_n-$hist_m.imputed.pval.signal.bigWig"
    done

    # return to the starting directory
    cd ..

done