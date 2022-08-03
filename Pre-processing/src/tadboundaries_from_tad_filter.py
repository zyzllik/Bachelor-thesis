# This script produces a file with TAD boundary windows and a file with TAD boundary window and neighbouring windows
# Windows exceeding chromosome boundaries are eliminated
# Gaps between TAD boundaries exceeding defined length (max_gap_length) are not considered
# Inputs: 
#   input_file_path
#   sample_name: e.g. cell_line
#   bin_number: OPTIONAL, Number of windows on each side from the TAD boundary window (e.g. bin_number=10 => 21 windows pro TAD boundary)
#   bin_size: OPTIONAL, Size of a window, only even numbers
#   Deafult: bin_number=10, bin_size=40000
# Output:
#   Number of detected boundaries that passed the filter

import sys
import json

# input_file_path = '/net/data.isilon/ag-cherrmann/projects/06_HIV_Microglia/data/tads/hg38/bed/K562_Rao_2014-raw_TADs.txt.bed'
# input_file_path = '/net/data.isilon/ag-cherrmann/echernova/IMR90/IMR90_TADs_lifted_hg38.bed'
input_file_path = sys.argv[1]
sample_name = sys.argv[2]
suffix = sys.argv[3]

if len(sys.argv)==6:
    bin_number = sys.argv[4]
    bin_size = sys.argv[5]
elif len(sys.argv)==4:
    bin_number = 10 
    bin_size = 40000
else:
    raise ValueError('Incorrect number of inputs')

output_file_path_windows = ('/net/data.isilon/ag-cherrmann/echernova/'
                            'tad_boundaries_{n_bins}_{bin_size}kb/'
                            '{sample}_TAD_boundaries_{n_bins}_{bin_size}kb{suffix}.txt.bed').format(n_bins=bin_number, 
                                                                                                     bin_size=bin_size//1000,
                                                                                                     sample=sample_name,
                                                                                                     suffix=suffix)
output_file_path_central = ('/net/data.isilon/ag-cherrmann/echernova/'
                            'tad_boundaries_{n_bins}_{bin_size}kb/'
                            '{sample}_TAD_boundaries_central_{bin_size}kb{suffix}.txt.bed').format(n_bins=bin_number, 
                                                                                                    bin_size=bin_size//1000,
                                                                                                    sample=sample_name,
                                                                                                    suffix=suffix)
 
max_gap_length = 4*10**6 # Maximal gap between TAD boundaries to be considered

input_f = open(input_file_path, 'r')
output_file_windows = open(output_file_path_windows, 'w')
output_file_central = open(output_file_path_central, 'w')

# chromosome_ends = {'chr1':248956422, 'chr2':242193529, 'chr3':198295559, 'chr4':190214555,
#     'chr5':181538259, 'chr6':170805979, 'chr7':159345973, 'chrX':156040895, 'chr8':145138636,
#     'chr9':138394717, 'chr11':135086622, 'chr10':133797422, 'chr12':133275309, 'chr13':114364328,
#     'chr14':107043718, 'chr15':101991189, 'chr16':90338345, 'chr17':83257441, 'chr18':80373285,
#     'chr20':64444167, 'chr19':58617616, 'chrY':57227415, 'chr22':50818468, 'chr21':46709983}

# with open('/net/data.isilon/ag-cherrmann/echernova/src/chrom_sizes.json', 'r') as chrom_sizes_file:
#     chromosome_ends = json.load(chrom_sizes_file)

if suffix=='_filtered':
    with open('/net/data.isilon/ag-cherrmann/echernova/src/chrom_sizes.json', 'r') as chrom_sizes_file:
        chromosome_ends = json.load(chrom_sizes_file)
elif suffix=='_imputed':
    chromosome_ends = {
        'chr1':249250621,
        'chr10':135534747,
        'chr11':135006516,
        'chr12':133851895,
        'chr13':115169878,
        'chr14':107349540,
        'chr15':102531392,
        'chr16':90354753,
        'chr17':81195210,
        'chr18':78077248,
        'chr19':59128983,
        'chr2':243199373,
        'chr20':63025520,
        'chr21':48129895,
        'chr22':51304566,
        'chr3':198022430,
        'chr4':191154276,
        'chr5':180915260,
        'chr6':171115067,
        'chr7':159138663,
        'chr8':146364022,
        'chr9':141213431,
        'chrX':155270560
    }
else:
    raise ValueError('Wrong suffix!')

line1 = input_f.readline()
counter = 0
for line2 in input_f:

    chr1 = line1.strip().split('\t')[0]
    chr2 = line2.strip().split('\t')[0]

    if chr1 == chr2:

        upstream_border = int(line1.strip().split('\t')[2])
        downstream_border = int(line2.strip().split('\t')[1])

        if (downstream_border-upstream_border) < max_gap_length:
            counter += 1
            border = (downstream_border+upstream_border)//2

            # Central window (bin)
            output_file_central.write("{chr}\t{u_border}\t{d_border}\n".format(
                chr =chr1, u_border=int(border-bin_size/2), d_border=int(border+bin_size/2)
            ))

            first_bin_u_border = border - bin_size//2 - bin_number*bin_size # upstream border of the first bin 
            for i in range(bin_number*2 + 1):
                bin_b = first_bin_u_border+bin_size*i
                if bin_b >= 0 and chromosome_ends[chr1]>=bin_b+bin_size:
                    # All windows (bins)
                    output_file_windows.write("{chr}\t{u_border}\t{d_border}\t{window}\t{tad_id}\n".format(chr=chr1, 
                                                                        u_border=int(bin_b), 
                                                                        d_border=int(bin_b+bin_size),
                                                                        window=i-bin_number,
                                        tad_id=counter))
    line1 = line2

print(counter) # Number of detected TAD boundaries that passed the filter 
input_f.close()
output_file_windows.close()
output_file_central.close()
