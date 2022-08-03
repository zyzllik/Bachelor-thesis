# -*- coding:utf-8 -*-
"""
Based on the work of Gan et al.
"""
from __future__ import print_function
from keras.models import Sequential
from keras.layers import Dropout, Dense, Activation, BatchNormalization
from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.metrics import roc_curve, auc
import numpy as np
from numpy.random import seed
import datetime
from pathlib import Path
import tensorflow as tf
from tensorflow.keras import layers

from src.utils import load_data

## -1 when excluded
#input_shape = load_data.img_cols * (load_data.img_rows -1)
batch_size = 128
num_classes = 2
epochs = 60


def mlp(input):
    """
    This function defines the MLP model

    Args:
        input: The input dataset.

    Returns:
        Compiled MLP model.
    """
    model = Sequential()
    model.add(Dense(512, input_shape=(input.shape[1],)))
    model.add(Activation('linear'))
    model.add(Dropout(0.6975))
    model.add(Dense(256))
    model.add(Activation('softplus'))
    model.add(Dropout(0.5153))
    model.add(BatchNormalization())
    model.add(Dense(512))
    model.add(Activation('linear'))
    model.add(Dropout(0.4252))
    model.add(BatchNormalization())
    model.add(Dense(1024))
    model.add(Activation('hard_sigmoid'))
    model.add(Dense(1))
    model.add(Activation('sigmoid'))

    model.compile(loss='binary_crossentropy', metrics=['accuracy'], optimizer='rmsprop')

    return model


def mlp_result(feature_data, result_folder, hist_list = False, exclude=False, date=None, input_type = 'csv'):
    """
    This function trains and tests MLP model.

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
    clf_mlp = KerasClassifier(build_fn=mlp, input=x_train_mlp, epochs=epochs, batch_size=batch_size, verbose=0)
    clf_mlp.fit(x_train_mlp, y_train_mlp, validation_data=(x_test_mlp, y_test_mlp))

    # Save the weights
    clf_mlp.save_weights(Path("/net/data.isilon/ag-cherrmann/echernova/lactuca/TAD-Lactuca/checkpoints/mlp_weights_{}".format(date)))
 
    # Test the model
    y_pred_mlp = clf_mlp.predict_proba(x_test_mlp)[:, 1]
    fpr_mlp, tpr_mlp, _ = roc_curve(y_test_mlp, y_pred_mlp)

    # Generate file names and save the results
    if hist_list is False:
        file_name = '{0}_MLP_full_model.npz'.format(date)
    else:
        file_name = '{0}_MLP_exclude_{1}'.format(date, '_'.join(hist_list))

    if not result_folder.exists():
        result_folder.mkdir()

    np.savez(result_folder/"{}".format(file_name), fpr = fpr_mlp, tpr = tpr_mlp, pred = y_pred_mlp, target = y_test_mlp)
    print('*******' * 3, '\n\t AUC {} = '.format(hist_list), auc(fpr_mlp, tpr_mlp), '\n', '*******' * 3)



if __name__ == '__main__':
    pass
