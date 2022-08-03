import pandas as pd
import numpy as np
from pathlib import Path
import itertools

def get_overlap_count(tads1, tads2, chrom, threshold=40000):
    # Calculates pairwise distance between center points of TAD boundaries
    # and compares it with a threshold. Returns number of overlaps
    center_arr1 = np.array(tads1[tads1['chr']==chrom]['center'])
    center_arr2 = np.array(tads2[tads2['chr']==chrom]['center'])
    distances = np.abs(np.subtract.outer(center_arr1, center_arr2))
    filtered_distances = distances < threshold
    counts1 = np.sum(filtered_distances, axis=1)
    counts2 = np.sum(filtered_distances, axis=0)
    return counts1, counts2

samples = [
    "A549", 
    "GM12878", 
    "HepG2", 
    "HMEC", 
    "HUVEC", 
    "IMR90", 
    "K562", 
    "NHEK",
    'H1-NPC',
    'Adrenal',
    'H1-TRO',
    'Ovary',
    'Aorta-STL002',
    'PANC1',
    'Bladder',
    'Pancreas',
    'Bowel-Small',
    'Hippocampus',
    'Psoas',
    'Caki2',
    'RPMI7951',
    'Cortex-DLPFC',
    'SJCRH30',
    'G401',
    'SKMEL5',
    'SKNDZ',
    'SKNMC',
    'GSE138767',
    'KBM7',
    'Spleen',
    'Liver-STL011',
    'T470',
    'H1-ESC',
    'LNCAP',
    'Ventricle-Left-STL003',
    'H1-MES',
    'Lung',
    'Ventricle-Right',
    'H1-MSC',
    'NCIH460'
]
suffix = '_imputed'
print(f"Working on data type: {suffix}")
folder_path = Path("/net/data.isilon/ag-cherrmann/echernova/tad_boundaries_10_40kb")
output_folder_path = Path("/net/data.isilon/ag-cherrmann/echernova/tad_bound_conservation_10_40kb")

all_tads = {}


# Load all TAD boundary files and prepare the data frames
print('Loading data')
for sample in samples:
    file_path = folder_path/f"{sample}_TAD_boundaries_10_40kb{suffix}.txt.bed"
    tads = pd.read_csv(file_path, names=['chr', 'start', 'end', 'window', 'tad_id'], sep='\t')
    tads = tads[tads['window']==0]
    tads['center'] = (tads['end']+tads['start'])//2
    tads.drop(columns=['window', 'start', 'end'], inplace=True)
    tads['conservation'] = 0
    all_tads[sample] = tads

# Get names of all chromosomes
chromosomes = np.unique(all_tads[samples[0]]['chr'])

# For each pair of cells calculate 
for (sample1, sample2) in itertools.combinations(samples, 2):
    print(f"Working on pair {sample1} and {sample2}")
    # Calculate distances for every chromosome separately
    for chrom in chromosomes:
        tads1 = all_tads[sample1]
        tads2 = all_tads[sample2]
        overlaps1, overlaps2 = get_overlap_count(tads1, tads2, chrom)
        # Add overlap count to the total overlap count
        tads1.loc[tads1['chr']==chrom, 'conservation'] = tads1[tads1['chr']==chrom]['conservation'] + overlaps1
        tads2.loc[tads2['chr']==chrom, 'conservation'] = tads2[tads2['chr']==chrom]['conservation'] + overlaps2
        # Save the updated dataframe in the dictionary
        all_tads[sample1] = tads1
        all_tads[sample2] = tads2

# Save data to csv
for sample in samples:
    all_tads[sample].to_csv(output_folder_path/f"{sample}_tad_bound_conservation{suffix}.csv")

