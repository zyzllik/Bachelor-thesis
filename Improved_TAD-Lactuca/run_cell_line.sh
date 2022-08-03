#!/bin/bash
# This file runs the MLP and RF model for all of the following cell lines
# Date: a word to be added to the result files (e.g. date)
date='0517'
for name in 'K562' 'HepG2' 'GM12878'
do
python src/tad_lactuca.py $name $date
done
