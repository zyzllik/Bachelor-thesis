"""
Contains functions for plotting the ROC and precision-recall curves.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn import metrics

def plot_roc(ax, file_path, data_name):
    """
    This function plots a ROC curve.

    Args:
        ax: axis to which to plot.
        file_path: Path to the .npz file including FPR, TPR and AUC.
        data_name: The label of the ROC curve 

    Returns:
        No return.

    """
    file = np.load(file_path, 'r')
    data_auc = metrics.auc(file['fpr'],file['tpr'])
    data_name = data_name + ' AUC = {0:.4f}'.format(data_auc)
    ax.plot(file['fpr'], file['tpr'], label = data_name)

def plot_roc_folder(folder, output_name, title):
    """
    This function plots a ROC curve for every .npz file in the folder
    to one figure and saves it.

    Args:
        folder: Path to the folder with .npz data.
        output_name: File to which to save the figure.
        title: Title of the figure.

    Returns:
        No return.

    """
    file_names = os.listdir(folder)
    fig, ax = plt.subplots()

    ax.set_prop_cycle(color=[
    '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a',
    '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
    '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
    '#17becf', '#9edae5'])

    for file in file_names:
        if file[:-4:-1] == 'zpn':
            path = folder / file
            if 'False' in file or 'full' in file:
                name = 'full model;'
            else:
                name = ', '.join(file.split('.')[0].split('_')[3:]) + ' excluded;'  
            plot_roc(ax, path, name)

    ax.set_xlabel('False positive rate')
    ax.set_ylabel('False negative rate')
    ax.set_title(title)
    ax.legend(prop={'size': 6})
    fig.tight_layout()
    fig.savefig(output_name)

def plot_precision_recall(ax, file_path, data_name):
    """
    This function plots a precision-recall curve.

    Args:
        ax: axis to which to plot.
        file_path: Path to the .npz file including FPR, TPR and AUPR.
        data_name: The label of the precision-recall curve 

    Returns:
        No return.

    """
    file = np.load(file_path, 'r')
    
    precision, recall, _ = metrics.precision_recall_curve(file['target'], file['pred'])
    aupr = metrics.average_precision_score(file['target'], file['pred'])
    data_name = data_name + ' AUPR = {0:.4f}'.format(aupr)
    ax.plot(recall, precision, label = data_name)

def plot_precision_recall_folder(folder, output_name, title):
    """
    This function plots a precision-recall curve for every .npz file in the folder
    to one figure and saves it.

    Args:
        folder: Path to the folder with .npz data.
        output_name: File to which to save the figure.
        title: Title of the figure.

    Returns:
        No return.

    """
    file_names = os.listdir(folder)
    fig, ax = plt.subplots()

    for file in file_names:
        if file[-4:] == '.npz':
            
            path = folder / file
            if 'False' in file or 'full' in file:
                name = 'full model;'
            else:
                name = '_'.join(file.split('-')[1:])[:-4] + ' excluded;'  
            plot_precision_recall(ax, path, name)

    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title(title)
    ax.legend(prop={'size': 6})
    fig.tight_layout()
    fig.savefig(output_name)


if __name__ == '__main__':
    pass
