#input_file_path = '/net/data.isilon/ag-cherrmann/projects/06_HIV_Microglia/data/tads/hg38/bed/K562_Rao_2014-raw_TADs.bed'
input_file_path = '/net/data.isilon/ag-cherrmann/echernova/IMR90/IMR90_TADs_lifted_hg38.bed'
output_file_path_windows = '/net/data.isilon/ag-cherrmann/echernova/tad_boundaries_10_40kb/IMR90_TAD_boundaries_10_40kb.txt.bed'
output_file_path_central = '/net/data.isilon/ag-cherrmann/echernova/tad_boundaries_10_40kb/IMR90_TAD_boundaries_central_40kb.txt.bed'
bin_number = 10
bin_size = 40000 # only even

input_f = open(input_file_path, 'r')
output_file_windows = open(output_file_path_windows, 'w')
output_file_central = open(output_file_path_central, 'w')

chromosome_ends = {'chr1':248956422, 'chr2':242193529, 'chr3':198295559, 'chr4':190214555,
    'chr5':181538259, 'chr6':170805979, 'chr7':159345973, 'chrX':156040895, 'chr8':145138636,
    'chr9':138394717, 'chr11':135086622, 'chr10':133797422, 'chr12':133275309, 'chr13':114364328,
    'chr14':107043718, 'chr15':101991189, 'chr16':90338345, 'chr17':83257441, 'chr18':80373285,
    'chr20':64444167, 'chr19':58617616, 'chrY':57227415, 'chr22':50818468, 'chr21':46709983}

line1 = input_f.readline()
counter = 0
for line2 in input_f:
    chr1 = line1.strip().split('\t')[0]
    chr2 = line2.strip().split('\t')[0]
    if chr1 == chr2:
        counter += 1
        upstream_border = int(line1.strip().split('\t')[2])
        downstream_border = int(line2.strip().split('\t')[1])
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

print("Processed TAD boundaries: {}".format(counter)) 
input_f.close()
output_file_windows.close()
output_file_central.close()
