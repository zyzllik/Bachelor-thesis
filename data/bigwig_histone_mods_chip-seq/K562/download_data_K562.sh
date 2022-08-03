#!/bin/bash
ids_list=(ENCFF525ZRM ENCFF405HIO ENCFF954JHK ENCFF834SEY ENCFF286WRJ ENCFF494WCA ENCFF849TDM ENCFF959YJV ENCFF901YVS ENCFF654SLZ ENCFF601JGK ENCFF605FAF)

mod_names=(H3K4me3 H3K27me3 H3K36me3 H3K4me1 H3K9ac H2AFZ H3K27ac H3K4me2 H3K79me2 H3K9me1 H3K9me3 H4K20me1)

for i in 0 1 2 3 4 5 6 7 8 9 10 11
do
	wget "https://www.encodeproject.org/files/${ids_list[i]}/@@download/${ids_list[i]}.bigWig"
	mv "${ids_list[i]}.bigWig" "HepG2_${mod_names[i]}_${ids_list[i]}"
done