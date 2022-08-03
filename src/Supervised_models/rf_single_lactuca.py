import os
import itertools
import json

import pandas as pd
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
import numpy as np
from numpy import random
import datetime
from pathlib import Path
import pickle
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier

def build_dataset(features, datapath):
    df = pd.read_csv(features).fillna(0)
    # first column: index from R, 2-4 columns: window coordinates --> need to be removed
    df.drop(['Unnamed: 0', 'chr', 'start', 'end'], axis=1, inplace=True)

    print('columns: {}'.format(df.columns))

    data = np.matrix(df.iloc[:, 4:])
    target = np.array(df.iloc[:, 3])
    print(data.shape, target.shape)
    # We do a 70-15-15 split
    x_train, x_rest, y_train, y_rest = train_test_split(data, target, train_size=0.7, random_state=0)
    x_valid, x_test, y_valid, y_test = train_test_split(x_rest, y_rest, train_size=0.5, random_state=0)
    print(f'Saving dataset to {datapath}!')
    np.savez(datapath, x_train=x_train,x_valid=x_valid, x_test=x_test,
                       y_train=y_train,y_valid=y_valid, y_test=y_test)
    print('Done')

def run_rf(datapath, dataset, params, seed=0, eval_model=True):
    
    # Load data
    data = np.load(dataset)
    train_x = data['x_train']
    train_y = data['y_train']
    valid_x = data['x_valid']
    valid_y = data['y_valid']

    # Train the model
    model = RandomForestClassifier(**params,
                                    # n_estimators=params['n_estimators'], 
                                    # max_depth=params['max_depth'],
                                    # min_samles_split = params['min_samples_split']
                                    random_state=seed)
    model.fit(train_x, train_y)

    # Save the model
    name = '_'.join(dataset.split('/')[-1].split('_')[:-2])
    with open(os.path.join(datapath,f"rf_weights_{name}_ne{params['n_estimators']}_md{params['max_depth']}_\
        ml{params['min_samples_leaf']}"), "wb") as out_file:
        pickle.dump(model, out_file)

    # Evaluate the model
    roc_auc = None
    if eval_model:
        roc_auc = model_performance(model, valid_x, valid_y, path=datapath, params=params, name=name)
    return roc_auc

def model_performance(model, x, y, path=None, params=None, name=None):

    # Test the model
    y_pred_rf = np.array(model.predict_proba(x))[:, 1]

    # One ROC curve pro class

    fpr_mlp, tpr_mlp, _ = roc_curve(y, y_pred_rf)
    roc_auc = auc(fpr_mlp, tpr_mlp)

    if path is not None:
        result = {'fpr': fpr_mlp, 'tpr': tpr_mlp, 'roc_auc': roc_auc}
        file_name = f"{name}_ne{params['n_estimators']}_md{params['max_depth']}_ml{params['min_samples_leaf']}_results.npz"
        file_path = os.path.join(path, file_name)
        np.savez(file_path, result)
        #plot_roc(fpr_mlp, tpr_mlp, roc_auc, name, params, path)
    print('*******' * 3, '\n\t AUC: ', roc_auc, '\n', '*******' * 3)
    return roc_auc

# def plot_roc(fpr, tpr, roc_auc, name, params, path):
#     fig, ax = plt.subplots()
#     classes = ['S1', 'S2', 'S3', 'not TAD']
#     for i in range(4):
#         ax.plot(fpr[i], tpr[i], label=f'{classes[i]}: AUC={roc_auc[i]:.3f}')
#     ax.set_title(f'ROC curve for {name}')
#     param_text = '\n'.join([f"{k}:{params[k]}" for k in params.keys()])
#     ax.text(0.7, 0.2, param_text)
#     ax.legend()
#     file_name = f"{name}_ne{params['n_estimators']}_b{params['bootstrap']}_md{params['max_depth']}_ml{params['min_samples_leaf']}_ms{params['min_samples_split']}_roc.png"
#     plt.savefig(os.path.join(path, file_name))

def main(features, params):
    name = features.split("/")[-1].split(".csv")[0]
    now = datetime.datetime.now().strftime("%Y%m%d")
    datapath = os.path.join(os.path.dirname(features), 'rf_single_allWindows_new')#{}'.format(now))
    dataset = os.path.join(datapath,f'{name}_dataset.npz')
    os.makedirs(datapath,exist_ok=True)
    if not os.path.exists(dataset):
        build_dataset(features, dataset)
    roc_auc = run_rf(datapath,
                    dataset,
                    params,
                    seed=3,
                    eval_model=True)
    return roc_auc

if __name__ == '__main__':
    print('------------------------------------')
    dir = '/net/data.isilon/ag-cherrmann/echernova/model_input_multi/'

    grid = {#'bootstrap': [True, False],
            'max_depth': [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, None],
            'min_samples_leaf': [1, 2, 4],
            #'min_samples_split': [2, 5, 10],
            'n_estimators': [200, 400, 600, 800, 1000]}

    param_permutations = [dict(zip(grid.keys(), v)) for v in itertools.product(*grid.values())]

    files = os.listdir(dir)

    auc_all=dict()

    for file in files:
        if file.endswith('allWindows.csv'):
            for parameters in param_permutations:
                print(f"File: {file}\n{parameters}")
                roc_auc = main(dir+file, params = parameters)

                all_auc_key = {'file': file, **parameters}
                auc_all[str(all_auc_key)] = roc_auc

    # Save all AUC scores
    name = files[0].split(".csv")[0]
    now = datetime.datetime.now().strftime("%Y%m%d")
    datapath = os.path.join(dir, 'rf_single_allWindows_new')
    print(datapath)
    os.makedirs(datapath,exist_ok=True)
    with open(os.path.join(datapath, 'all_auc.json'), "w") as output:
        json.dump(auc_all, output)
    
    