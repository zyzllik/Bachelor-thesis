from pathlib import Path 
import json

window_size = 40000
output_path = Path('/net/data.isilon/ag-cherrmann/echernova/integrative_nmf/220701_all_windows_filtered.txt.bed')

with open('/net/data.isilon/ag-cherrmann/echernova/src/chrom_sizes.json', 'r') as chrom_sizes_file:
    chromosome_ends = json.load(chrom_sizes_file)

output_file = open(output_path, 'w')

for chrom in chromosome_ends.keys():
    chrom_end = chromosome_ends[chrom]
    chrom_start = 0
    while chrom_start + window_size < chrom_end-1:
        output_file.write(f'{chrom}\t{chrom_start}\t{chrom_start+window_size}\n')
        chrom_start += window_size
    output_file.write(f'{chrom}\t{chrom_start}\t{chrom_end}\n')

output_file.close()

