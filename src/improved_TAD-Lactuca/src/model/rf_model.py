# -*- coding:utf-8 -*-
"""
Based on the work of Gan et al.
"""
from __future__ import print_function
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import RandomForestClassifier

from src.utils import load_data

import numpy as np
from pathlib import Path
import pickle


def rf_result(feature_data, result_folder, hist_list = False, 
              exclude=False, date=None, input_type = 'csv'):
    """
    This function trains and tests RF model

    Args:
        feature_data: Tuple including path to the input data of class 1 and 0.
        result_folder: Path to the folder where results should be saved.
        hist_list: List of histone modifications that should be excluded 
                   or included based on the value of exclude. If False all
                   features are included.
        exclude: If True features in the hist_list are excluded from the dataset, 
                 if False features in the hist_list are included.
        date: Date or other word to be added to the file names.
        input_type: Type of the input data (csv or xlsx).

    Returns:
        Nothing.

    """
    
    # Load data
    if input_type=='csv':
        feature_y, feature_n = feature_data
        (x_train_mlp, y_train_mlp), (x_test_mlp, y_test_mlp) = load_data.csv_load(feature_y, feature_n, hist_mod_list=hist_list, exclusion=exclude)
    elif input_type=='xlsx':
        (x_train_mlp, y_train_mlp), (x_test_mlp, y_test_mlp) = load_data.xlsx_load(feature_data, hist_mod_list=hist_list, exclusion=exclude)
    else:
        print("Wrong input type!")

    # Train the model
    clf_rf = RandomForestClassifier(n_estimators=500)
    clf_rf.fit(x_train_mlp, y_train_mlp)

    # Save the model
    pickle.dump(clf_rf, 
                open(Path("/net/data.isilon/ag-cherrmann/echernova/lactuca/TAD-Lactuca/checkpoints/rf_weights_{}".format(date)), "wb"))

    # Test the model
    y_pred_rf = clf_rf.predict_proba(x_test_mlp)[:, 1]
    fpr_rf, tpr_rf, _ = roc_curve(y_test_mlp, y_pred_rf)

    # Generate file names and save the results
    if hist_list is False:
        file_name = '{0}_RF_full_model.npz'.format(date)
    else:
        file_name = '{0}_RF_exclude_{1}'.format(date, '_'.join(hist_list))

    if not result_folder.exists():
        result_folder.mkdir()
    
    np.savez(result_folder/"{}".format(file_name), fpr = fpr_rf, tpr = tpr_rf, pred = y_pred_rf, target = y_test_mlp)
    print('*******' * 3, '\n\t AUC = ', auc(fpr_rf, tpr_rf), '\n', '*******' * 3)

if __name__ == '__main__':
    pass
