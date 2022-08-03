#!/bin/bash
url_list=(https://www.encodeproject.org/files/ENCFF003DXG/@@download/ENCFF003DXG.bigWig https://www.encodeproject.org/files/ENCFF919DOR/@@download/ENCFF919DOR.bigWig https://www.encodeproject.org/files/ENCFF831WJT/@@download/ENCFF831WJT.bigWig https://www.encodeproject.org/files/ENCFF935EGN/@@download/ENCFF935EGN.bigWig https://www.encodeproject.org/files/ENCFF340JIF/@@download/ENCFF340JIF.bigWig https://www.encodeproject.org/files/ENCFF564KBE/@@download/ENCFF564KBE.bigWig https://www.encodeproject.org/files/ENCFF627OKN/@@download/ENCFF627OKN.bigWig https://www.encodeproject.org/files/ENCFF931USZ/@@download/ENCFF931USZ.bigWig https://www.encodeproject.org/files/ENCFF599TRR/@@download/ENCFF599TRR.bigWig https://www.encodeproject.org/files/ENCFF683HCZ/@@download/ENCFF683HCZ.bigWig https://www.encodeproject.org/files/ENCFF479XIQ/@@download/ENCFF479XIQ.bigWig)

ids_list=(ENCFF003DXG ENCFF919DOR ENCFF831WJT ENCFF935EGN ENCFF340JIF ENCFF564KBE ENCFF627OKN ENCFF931USZ ENCFF599TRR ENCFF683HCZ ENCFF479XIQ)

mod_names=(H3K4me3 H3K27me3 H3K36me3 H2AFZ H3K27ac H3K4me1 H3K4me2 H3K79me2 H3K9ac H3K9me3 H4K20me1)

for i in 0 1 2 3 4 5 6 7 8 9 10 11
do
	wget ${url_list[i]}
	mv "${ids_list[i]}.bigWig" "GM12878_${mod_names[i]}_${ids_list[i]}.bigWig"

done
