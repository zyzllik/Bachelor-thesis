#!/bin/bash
ids_list=(ENCFF105IIK ENCFF425LVX ENCFF074PND ENCFF431PAM ENCFF121EJY)

mod_names=(H3K27ac H3K9me3 H3K4me3 H3K4me1 H3K4me2)

for i in 0 1 2 3 4
do
	wget "https://www.encodeproject.org/files/${ids_list[i]}/@@download/${ids_list[i]}.bigWig"
	mv "${ids_list[i]}.bigWig" "HepG2_${mod_names[i]}_${ids_list[i]}"
done