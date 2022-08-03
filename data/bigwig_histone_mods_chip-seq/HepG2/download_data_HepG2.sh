#!/bin/bash
url_list=(https://www.encodeproject.org/files/ENCFF736LHE/@@download/ENCFF736LHE.bigWig https://www.encodeproject.org/files/ENCFF598TWA/@@download/ENCFF598TWA.bigWig https://www.encodeproject.org/files/ENCFF550PXE/@@download/ENCFF550PXE.bigWig https://www.encodeproject.org/files/ENCFF253PND/@@download/ENCFF253PND.bigWig https://www.encodeproject.org/files/ENCFF795ONN/@@download/ENCFF795ONN.bigWig https://www.encodeproject.org/files/ENCFF576YVM/@@download/ENCFF576YVM.bigWig https://www.encodeproject.org/files/ENCFF057BKO/@@download/ENCFF057BKO.bigWig https://www.encodeproject.org/files/ENCFF655XBP/@@download/ENCFF655XBP.bigWig https://www.encodeproject.org/files/ENCFF040RHK/@@download/ENCFF040RHK.bigWig https://www.encodeproject.org/files/ENCFF754ROM/@@download/ENCFF754ROM.bigWig https://www.encodeproject.org/files/ENCFF330AIV/@@download/ENCFF330AIV.bigWig)

ids_list=(ENCFF736LHE ENCFF598TWA ENCFF550PXE ENCFF253PND ENCFF795ONN ENCFF576YVM ENCFF057BKO ENCFF655XBP ENCFF040RHK ENCFF754ROM ENCFF330AIV)

mod_names=(H3K4me3 H3K27me3 H3K36me3 H2AFZ H3K27ac H3K4me1 H3K4me2 H3K79me2 H3K9ac H3K9me3 H4K20me1)

for i in 0 1 2 3 4 5 6 7 8 9 10 11
do
	wget ${url_list[i]}
	mv "${ids_list[i]}.bigWig" "HepG2_${mod_names[i]}_${ids_list[i]}"
done