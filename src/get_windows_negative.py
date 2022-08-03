import sys

sample = sys.argv[1]
suffix = sys.argv[2]

input_file_path = ('/net/data.isilon/ag-cherrmann/echernova/'
                    'tad_boundaries_10_40kb/'
                    '{sample}_negatives_central_40kb{suffix}.txt.bed').format(sample=sample, suffix=suffix)
output_file_path = ('/net/data.isilon/ag-cherrmann/echernova/'
                    'tad_boundaries_10_40kb/'
                    '{sample}_negatives_10_40kb{suffix}.txt.bed').format(sample=sample, suffix=suffix)

bin_number = 10
bin_size = 40000 # only even

input_f = open(input_file_path, 'r')
output_file = open(output_file_path, 'w')

# chromosome_ends = {'chr1':248956422, 'chr2':242193529, 'chr3':198295559, 'chr4':190214555,
#     'chr5':181538259, 'chr6':170805979, 'chr7':159345973, 'chrX':156040895, 'chr8':145138636,
#     'chr9':138394717, 'chr11':135086622, 'chr10':133797422, 'chr12':133275309, 'chr13':114364328,
#     'chr14':107043718, 'chr15':101991189, 'chr16':90338345, 'chr17':83257441, 'chr18':80373285,
#     'chr20':64444167, 'chr19':58617616, 'chrY':57227415, 'chr22':50818468, 'chr21':46709983}

# Imputed chrom ends
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

counter = 0
for line in input_f:
    counter += 1
    chrom = line.strip().split('\t')[0]
    central_u_border =int(line.strip().split('\t')[1])
    central_d_border = int(line.strip().split('\t')[2])

    border = (central_u_border+central_d_border)//2
    first_bin_u_border = border - bin_size//2 - bin_number*bin_size # upstream border of the first bin 

    for i in range(bin_number*2+1):
        bin_b = first_bin_u_border+bin_size*i
        if bin_b >= 0 and chromosome_ends[chrom]>=bin_b+bin_size:
            output_file.write("{chr}\t{u_border}\t{d_border}\t{window}\t{tad_id}\n".format(chr=chrom, 
                                                                u_border=bin_b, 
                                                                d_border=bin_b+bin_size,
                                                                window=i-bin_number,
								                                tad_id=counter))
print("Processed TAD boundaries: {}".format(counter)) 
input_f.close()
output_file.close()
