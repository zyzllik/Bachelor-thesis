#!/bin/sh

source ~/miniconda3/bin/activate
conda activate testenv
python /net/data.isilon/ag-cherrmann/echernova/src/tadboundaries_from_tad_filter.py
