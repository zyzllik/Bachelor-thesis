#!/bin/bash

directory="/net/data.isilon/ag-cherrmann/projects/06_HIV_Microglia/data/tads/hg38/bed/"
# tad_input_files=(
#     'A549_raw-merged_TADs.txt.bed'
#     'GM12878_Rao_2014-raw_TADs.txt.bed'
#     'HepG2.ENCODE3-raw_modified.bed'
#     'HMEC_Rao_2014-raw_TADs.txt.bed'
#     'HUVEC_Rao_2014-raw_TADs.txt.bed'
#     'IMR90_Rao_2014-raw_TADs.txt.bed'
#     'K562_Rao_2014-raw_TADs.txt.bed'
#     'NHEK_Rao_2014-raw_TADs.txt.bed'
# )

# samples=(
#     "A549"
#     "GM12878"
#     "HepG2"
#     "HMEC"
#     "HUVEC"
#     "IMR90"
#     "K562"
#     "NHEK"
# )
tad_input_files=(
    'H1-NPC_Dixon_2015-raw_TADs.txt.bed'
    'Adrenal_Schmitt2016-raw_TADs.txt.bed'
    'H1-TRO_Dixon_2015-raw_TADs.txt.bed'
    'Ovary_Schmitt2016-raw_TADs.txt.bed'
    'Aorta_STL002_Leung_2015-raw_TADs.txt.bed'
    'PANC1_raw-merged_TADs.txt.bed'
    'Bladder_Schmitt2016-raw_TADs.txt.bed'
    'Pancreas_Schmitt2016-raw_TADs.bed'
    'Bowel_Small_Schmitt2016-raw_TADs.txt.bed'
    'Hippocampus_Schmitt2016-raw_TADs.bed'
    'Psoas_Schmitt2016-raw_TADs.txt.bed'
    'Caki2_raw-merged_TADs.txt.bed'
    'RPMI7951_raw-merged_TADs.txt.bed'
    'Cortex_DLPFC_Schmitt2016-raw_TADs.bed'
    'SJCRH30_raw-merged_TADs.txt.bed'
    'G401_raw-merged_TADs.txt.bed'
    'SKMEL5_raw-merged_TADs.txt.bed'
    'SKNDZ_raw-merged_TADs.txt.bed'
    'SKNMC_raw-merged_TADs.txt.bed'
    'GSE138767_CD4_merged_TADS.bed'
    'KBM7_Rao_2014-raw_TADs.txt.bed'
    'Spleen_Schmitt2016-raw_TADs.txt.bed'
    'Liver_STL011_Leung_2015-raw_TADs.txt.bed'
    'T470_raw-merged_TADs.txt.bed'
    'H1-ESC_Dixon_2015-raw_TADs.txt.bed'
    'LNCAP_raw-merged_TADs.txt.bed'
    'VentricleLeft_STL003_Leung_2015-raw_TADs.txt.bed'
    'H1-MES_Dixon_2015-raw_TADs.txt.bed'
    'Lung_Schmitt2016-raw_TADs.txt.bed'
    'Ventricle_Right_Schmitt2016-raw_TADs.txt.bed'
    'H1-MSC_Dixon_2015-raw_TADs.txt.bed'
    'NCIH460_raw-merged_TADs.txt.bed'
)

samples=(
    'H1-NPC'
    'Adrenal'
    'H1-TRO'
    'Ovary'
    'Aorta-STL002'
    'PANC1'
    'Bladder'
    'Pancreas'
    'Bowel-Small'
    'Hippocampus'
    'Psoas'
    'Caki2'
    'RPMI7951'
    'Cortex-DLPFC'
    'SJCRH30'
    'G401'
    'SKMEL5'
    'SKNDZ'
    'SKNMC'
    'GSE138767'
    'KBM7'
    'Spleen'
    'Liver-STL011'
    'T470'
    'H1-ESC'
    'LNCAP'
    'Ventricle-Left-STL003'
    'H1-MES'
    'Lung'
    'Ventricle-Right'
    'H1-MSC'
    'NCIH460'
)

suffix='_filtered'


# !!!! 01 Changed (imputed bigwig files)

# Activate the right conda enviroment
source ~/miniconda3/bin/activate
conda activate testenv

n_samples=${#samples[@]}

for ((i=0; i<$n_samples; i++))
do
    sample=${samples[$i]}
    tad_file="$directory${tad_input_files[$i]}"

    echo "------------ Processing sample: $sample ------------"

    echo "## Positives ##"

    ### Get TAD boundaries from the TADs ###
    # The python script creates a file with the central window of the TAD boundary
    # and a file with central and neighbouring windows of the TAD boundary
    # Number of TAD boundaries is saved in n_tads
    echo "Extract TAD boundaries from the TAD file"
    echo $tad_file
    n_tads=`python /net/data.isilon/ag-cherrmann/echernova/src/tadboundaries_from_tad_filter.py $tad_file $sample $suffix`
    echo "Number of TAD boundaries detected: $n_tads"

    # ### Get binned features in the defined windows ###
    # echo "Binning features"
    # path_TAD_boundaries="/net/data.isilon/ag-cherrmann/echernova/tad_boundaries_10_40kb/${sample}_TAD_boundaries_10_40kb${suffix}.txt.bed"
    # sh /net/data.isilon/ag-cherrmann/echernova/src/01-definingwindows.sh $sample $path_TAD_boundaries $suffix
    
    # echo "## Negatives ##"

    # ### Generate central windows for negative samples ###
    # echo "Generating synthetic negative central TAD boundary windows"
    # sh /net/data.isilon/ag-cherrmann/echernova/src/n01_get_central.sh $n_tads $sample $suffix

    # ### Get neighbouring windows for the negative samples ###
    # echo "Calculating neigbouring windows"
    # python /net/data.isilon/ag-cherrmann/echernova/src/get_windows_negative.py $sample $suffix

    # ### Get binned features in the negative windows
    # echo "Binning features"
    # path_to_neg="/net/data.isilon/ag-cherrmann/echernova/tad_boundaries_10_40kb/${sample}_negatives_10_40kb${suffix}.txt.bed"
    # sh /net/data.isilon/ag-cherrmann/echernova/src/01-definingwindows.sh $sample $path_to_neg "${suffix}_neg"

done