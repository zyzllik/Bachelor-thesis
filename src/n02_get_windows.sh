#!/bin/sh

source ~/miniconda3/bin/activate
conda activate testenv
python /net/data.isilon/ag-cherrmann/echernova/src/get_windows_negative.py
