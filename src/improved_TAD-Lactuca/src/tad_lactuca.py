# -*- coding:utf-8 -*-
"""
Based on the work of Gan et al.
Improved to compare the full and the reduced models. 
Expects the name of the cell line and a further string (for example the date) to be passed on.
"""
import sys
from pathlib import Path

workspace = './'
sys.path.append(workspace)


from model import mlp_model, rf_model
from utils.plots import *

if __name__ == '__main__':


    # Folder containig the input data
    input_folder = './'
    # Folder where the output of the model (ROC curve) should be saved
    output_folder = './'
    # Testing the full model and the model with features available for microglia
    exclude = False # the list is included for each list in data_list
    data_list = [False, ['H3K9me3', 'H3K4me1', 'H3K36me3', 'H3K27me3', 'H3K27ac', 'CTCF']] 
        # full (False) and reduced model (the features from the second list) will be trained



    # Load inputs
    if len(sys.argv) == 3:
        cell_line = sys.argv[1]
        date = sys.argv[2]
        feature_pos= input_folder/"{0}/positives_{1}_filtered.csv".format(cell_line, cell_line)
        feature_neg= input_folder/"/{0}/negatives_{1}_filtered.csv".format(cell_line, cell_line)
    else:
        print("Incorrect number of inputs!")
    
    # MLP model
    print("MLP...")
    result_folder_mlp = output_folder/'/{0}'.format(cell_line)
    for exluded_features in data_list:
        print(exluded_features)
        mlp_model.mlp_result((feature_pos, feature_neg), result_folder_mlp, hist_list=exluded_features, exclude=exclude, date=date, input_type='csv')
    plot_roc_folder(result_folder_mlp, result_folder_mlp/'{0}_mlp_ROC_curve_{1}_only_available_mods.png'.format(date, cell_line), 'ROC comparison: MLP on {0}'.format(cell_line))
    
    # RF model
    print("RF...")
    result_folder_rf = output_folder/'/{0}'.format(cell_line)
    for exluded_features in data_list:
        print(exluded_features)
        rf_model.rf_result((feature_pos, feature_neg), result_folder_rf, hist_list=exluded_features, exclude=exclude, date=date, input_type='csv')
    plot_roc_folder(result_folder_rf, result_folder_rf/'{0}_rf_ROC_curve_{1}_only_available_mods.png'.format(date, cell_line), 'ROC comparison: RF on {0}'.format(cell_line))
