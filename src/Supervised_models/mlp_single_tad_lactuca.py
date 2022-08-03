import os
import json

import pandas as pd
from sklearn.metrics import roc_curve, auc
import numpy as np
from numpy import random
import datetime
from pathlib import Path

import tensorflow as tf
from sklearn.model_selection import train_test_split
from tensorflow.keras import layers

import matplotlib.pyplot as plt


# batch_size = 128
# num_classes = 2
# epochs = 60

def build_dataset(features, datapath):
    df = pd.read_csv(features).fillna(0)
    # first column: index from R, 2-4 columns: window coordinates --> need to be removed
    df.drop(['Unnamed: 0', 'chr', 'start', 'end'], axis=1, inplace=True)

    print('columns: {}'.format(df.columns))

    data = np.matrix(df.iloc[:, 4:])

    target = np.array(df.iloc[:, 3]) # Only the not bound column in this case
    # We do a 70-15-15 split
    x_train, x_rest, y_train, y_rest = train_test_split(data, target, train_size=0.7, random_state=0)
    x_valid, x_test, y_valid, y_test = train_test_split(x_rest, y_rest, train_size=0.5, random_state=0)
    print(f'Saving dataset to {datapath}!')
    np.savez(datapath, x_train=x_train,x_valid=x_valid, x_test=x_test,
                       y_train=y_train,y_valid=y_valid, y_test=y_test)
    print('Done')


def build_mlp(input,layersize=512,dropout=0.0,lr=0.001):

    model = tf.keras.models.Sequential([
        layers.Dense(layersize, input_shape=(input.shape[1],), activation='relu'),
        layers.Dropout(dropout),
        layers.Dense(layersize, activation='relu'),
        layers.Dropout(dropout),
        layers.BatchNormalization(), 
        layers.Dense(layersize, activation='relu'),
        layers.Dropout(dropout),
        layers.BatchNormalization(),
        layers.Dense(layersize, activation='relu'),
        layers.Dense(1, activation='sigmoid')
    ])
    print(lr)
    model.compile(loss='binary_crossentropy', metrics=['accuracy'], optimizer=tf.keras.optimizers.Adam(lr=lr))  # optimizer=adam
    return model


def run_mlp(datapath, dataset, seed=0,logging_dir='logs/model',
            batch_size=128,epochs=10,learning_rate=0.001,
            dropout=0.0, layersize=512, eval_model=True):
    random.seed(seed)

    # Load data
    data = np.load(dataset)
    train_x = data['x_train']
    train_y = data['y_train']
    valid_x = data['x_valid']
    valid_y = data['y_valid']

    # Train the model
    name = '_'.join(dataset.split('/')[-1].split('_')[:-2])
    #now = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    log_dir = f'{logging_dir}_{name}_lr{learning_rate}_drop{dropout}_ls{layersize}'
    print("Saving logs to: ", log_dir)
    tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1)

    model = build_mlp(train_x,lr=learning_rate, dropout=dropout, layersize=layersize)
    model.fit(x=train_x,
              y=train_y,
              epochs=epochs,
              batch_size=batch_size,
              validation_data=(valid_x, valid_y),
              callbacks=[tensorboard_callback])

    # Save the weights
    model.save_weights(os.path.join(datapath,f"weights/mlp_weights_{name}_ls{learning_rate}_drop{dropout}_ls{layersize}"))

    if eval_model:
        model_performance(model,valid_x,valid_y,datapath, name, learning_rate, dropout=dropout, layersize=layersize)




def model_performance(model,x,y,path=None, name=None, lr=None, dropout=None, layersize=None):

    # Test the model
    y_pred_mlp = model.predict(x)[:,0]
    print(y_pred_mlp)

    # One ROC curve pro class
    fpr_mlp = dict()
    tpr_mlp = dict()
    roc_auc = dict()
    fpr_mlp, tpr_mlp, _ = roc_curve(y, y_pred_mlp)
    roc_auc = auc(fpr_mlp, tpr_mlp)

    if path is not None:
        result = {'fpr': fpr_mlp, 'tpr': tpr_mlp, 'roc_auc': roc_auc}
        file_path = os.path.join(path, f'{name}_lr{lr}_drop{dropout}_ls{layersize}_results.npz')
        np.savez(file_path, result)
        #np.savez(os.path.join(path,'results.npz'), fpr=fpr_mlp, tpr=tpr_mlp, pred=y_pred_mlp, target=y)
        #plot_roc(fpr_mlp, tpr_mlp, roc_auc, name, lr, dropout, layersize, path)
    print('*******' * 3, '\n\t AUC: ', roc_auc, '\n', '*******' * 3)


def plot_roc(fpr, tpr, roc_auc, name, lr, dropout, layersize, path):
    fig, ax = plt.subplots()
    classes = ['S1', 'S2', 'S3', 'not TAD']
    # for i in range(4):
    #     ax.plot(fpr, tpr, label=f'{classes[i]}: AUC={roc_auc[i]:.3f}')
    cell_name = name.split('_')[0]
    quantile = name.split('_')[1]
    ax.set_title(f'ROC curve for {cell_name}\nQuantile: {quantile[1:]}, learning rate: {lr},\ndropout: {dropout}, layersize: {layersize}')
    ax.legend()
    plt.savefig(os.path.join(path, f'{name}_lr{lr}_drop{dropout}_ls{layersize}_roc.png'))



def main(features, lr=0.001, dropout=0.0, layersize=512):
    name = features.split("/")[-1].split(".csv")[0]
    now = datetime.datetime.now().strftime("%Y%m%d")
    datapath = os.path.join(os.path.dirname(features), 'mlp_single_allWindows')#{}'.format(now))
    dataset = os.path.join(datapath,f'{name}_dataset.npz')
    os.makedirs(datapath,exist_ok=True)
    if not os.path.exists(dataset):
        build_dataset(features, dataset)
    n_epochs = 100
    run_mlp(datapath,
            dataset,
            epochs=n_epochs,
            learning_rate=lr,
            dropout=dropout,
            layersize=layersize,
            eval_model=True,
            logging_dir=os.path.join(datapath,'logs/model'))



if __name__ == '__main__':
    print('------------------------------------')
    dir = '/net/data.isilon/ag-cherrmann/echernova/model_input_multi/'
    learning_rates = [0.0001, 0.0005, 0.001, 0.005, 0.01]
    dropout_list = [0.0, 0.1, 0.2, 0.3, 0.4]
    layersize_list = [128, 256, 512]
    files = os.listdir(dir)
    for file in files:
        if file.endswith('allWindows.csv'):
            for lr in learning_rates:
                for  dropout in dropout_list:
                    for layersize in layersize_list:
                        print(f"File: {file}\tLearning rate: {lr}\tDropout: {dropout}\tLayersize: {layersize}")
                        main(dir+file, lr=lr, dropout=dropout, layersize=layersize)