# -*- coding:utf-8 -*-
"""
Based on work by Gan et al.
Improved functions for data loading from .csv and .xlsx files.
"""

from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np


def csv_load(feature_y, feature_n, hist_mod_list=False, exclusion=False):
    """
    This function loads data from the .csv file, excludes some columns if needed 
    and returns dataset split in train (70%) and test (30%).

    Args:
        feature_y: Path to the data of class 1.
        feature_n: Path to the data of class 0.
        hist_mod_list: List of features to be excluded or included
                       based on the value of exclusion. If False
                       all features will be included
        exclusion: If True the features in the hist_mod_list will be
                   excluded, if False the features in the hist_mod_list
                   will be included.

    Returns:
        (x_train, y_train), (x_test, y_test): The dataset split in train and test.

    """

    df_y = pd.read_csv(feature_y).fillna(0)
    df_n = pd.read_csv(feature_n).fillna(0)
    # first column: index from R, second column:tad_ids --> need to be removed
    df_y.drop(['Unnamed: 0', 'tad_id'], axis=1,  inplace = True)
    df_n.drop(['Unnamed: 0', 'tad_id'], axis=1,  inplace = True)
    # append label column
    df_y['label'] = 1
    df_n['label'] = 0
    df = pd.concat([df_y, df_n], axis=0)

    if hist_mod_list is not False:
        if exclusion:
            drop_columns = []
            for hist_mod in hist_mod_list:
                drop_columns += [col_name for col_name in list(df.columns) if hist_mod in col_name]
            df = df.drop(drop_columns, axis=1)
        elif not exclusion: # aka inclusion
            keep_columns = ['chr', 'start', 'end']
            for hist_mod in hist_mod_list:
                keep_columns += [col_name for col_name in list(df.columns) if hist_mod in col_name]
            keep_columns += ['label']
            df = df[keep_columns]
        
    print('columns: {}'.format(df.columns))
    
    data = np.matrix(df.iloc[:,3:-1])
    target = np.array(df['label'])

    x_train, x_test, y_train, y_test = train_test_split(data, target, train_size=0.7, random_state=49)

    return (x_train, y_train), (x_test, y_test)

def xlsx_load(feature, hist_mod_list=False, exclusion=False):
    """
    This function loads data from the .xlsx file, excludes some columns if needed 
    and returns dataset split in train (70%) and test (30%).

    Args:
        feature_y: Path to the data of class 1.
        feature_n: Path to the data of class 0.
        hist_mod_list: List of features to be excluded or included
                       based on the value of exclusion. If False
                       all features will be included
        exclusion: If True the features in the hist_mod_list will be
                   excluded, if False the features in the hist_mod_list
                   will be included.

    Returns:
        (x_train, y_train), (x_test, y_test): The dataset split in train and test.

    """

    df_y = pd.read_excel(feature, sheet_name='y').fillna(0)
    df_n = pd.read_excel(feature, sheet_name='n').fillna(0)
    df = df_y.append(df_n)

    if hist_mod_list is not False:
        if exclusion:
            drop_columns = []
            for hist_mod in hist_mod_list:
                drop_columns += [col_name for col_name in list(df.columns) if hist_mod in col_name]
            df = df.drop(drop_columns, axis=1)
        elif not exclusion: # aka inclusion
            keep_columns = ['chr', 'start', 'end']
            for hist_mod in hist_mod_list:
                keep_columns += [col_name for col_name in list(df.columns) if hist_mod in col_name]
            keep_columns += ['label']
            df = df[keep_columns]
        
    print('columns: {}'.format(df.columns))
    
    index = [i for i in range(4, df.shape[1])]
    data = np.matrix(df.iloc[:, index])
    target = np.array(df.iloc[:, 3])
    x_train, x_test, y_train, y_test = train_test_split(data, target, train_size=0.7, random_state=49)
    return (x_train, y_train), (x_test, y_test)



if __name__ == '__main__':
    pass
