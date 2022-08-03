#!/bin/bash
#url_list=(https://www.encodeproject.org/files/ENCFF221MJG/@@download/ENCFF221MJG.bigWig https://www.encodeproject.org/files/ENCFF993MGL/@@download/ENCFF993MGL.bigWig https://www.encodeproject.org/files/ENCFF376ZIM/@@download/ENCFF376ZIM.bigWig https://www.encodeproject.org/files/ENCFF355BSW/@@download/ENCFF355BSW.bigWig https://www.encodeproject.org/files/ENCFF066MEE/@@download/ENCFF066MEE.bigWig https://www.encodeproject.org/files/ENCFF699OAR/@@download/ENCFF699OAR.bigWig https://www.encodeproject.org/files/ENCFF366BVS/@@download/ENCFF366BVS.bigWig https://www.encodeproject.org/files/ENCFF519FHW/@@download/ENCFF519FHW.bigWig https://www.encodeproject.org/files/ENCFF105FHL/@@download/ENCFF105FHL.bigWig)

#ids_list=(ENCFF221MJG ENCFF993MGL ENCFF376ZIM ENCFF355BSW ENCFF066MEE ENCFF699OAR ENCFF366BVS ENCFF519FHW ENCFF105FHL)

#mod_names=(H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K9me3 H3K27ac H3K27me3 H3K36me3 CTCF)

url_list=(
	https://www.encodeproject.org/files/ENCFF859KTR/@@download/ENCFF859KTR.bigWig
	https://www.encodeproject.org/files/ENCFF295EWC/@@download/ENCFF295EWC.bigWig
	https://www.encodeproject.org/files/ENCFF050ZTH/@@download/ENCFF050ZTH.bigWig
)

ids_list=(ENCFF859KTR ENCFF295EWC ENCFF050ZTH)

mod_names=(H2AFZ H4K20me1 H3K79me2)

for i in 0 1 2 3 4 5 6 7 8 9 10 11
do
	wget --no-check-certificate ${url_list[i]}
	mv "${ids_list[i]}.bigWig" "GM12878_${mod_names[i]}_${ids_list[i]}.bigWig"

done
