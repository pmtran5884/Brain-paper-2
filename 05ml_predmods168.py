# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 10:25:15 2019

@author: ptran
"""

import os
os.chdir('/Users/ptran/Box Sync/TCGA Cancer Classification Project/Brain ns processing 20190814')
from loading import get_optimal_genes, clean_name
import numpy as np
import pandas as pd
from datetime import datetime
from sklearn.svm import SVC
from sklearn.metrics import classification_report, silhouette_score
from sklearn.model_selection import train_test_split, cross_val_score
from random import shuffle
from scipy.stats import mode

# hyperparameters
num_models = 1000
min_genes_3grps = 25
num_gene_pool_3grps = 100

# load data
data_file = '03bBrain_MCGandTCGA_CombinedData_20190814_sampsrows.csv'
data = pd.read_csv(data_file, index_col=[0])
print(data.shape)
print(data.columns)
print(data.head())

# load genes and update data
genes_hall_3grps = pd.read_csv('She_selected_168.csv').squeeze()
genes_common_3grps = set(data.columns).intersection(genes_hall_3grps["Gene"])
data_3grps = data.reindex(genes_common_3grps, axis=1)

# assign targets
clin = pd.read_csv('03bBrain_MCGandTCGA_PhenoData_20190821.csv',
                   index_col=[0])
target_3grps = clin.DBU_4grp

print(target_3grps.value_counts())

# remove ambi and mcg
target_3grps = target_3grps.loc[clin.Study=='TCGA']
target_3grps = target_3grps.loc[target_3grps != 'Ambi']
data_3grps = data_3grps.reindex(target_3grps.index)

# clean data
print(data_3grps.isnull().any().any())
print(data_3grps.shape)

# generate models for 3grps
now = clean_name(str(datetime.now()))
iteration_folder = 'validation_hallmarks_{}'.format(now)
os.mkdir(iteration_folder)
os.chdir(iteration_folder)
os.mkdir('3grps')
#os.mkdir('sub2')
os.mkdir('dropout_shuffled_models')
os.mkdir('dropout_good_model_umaps')

#playtime

random_genes_dict = {}
for model in range(1, num_models + 1):
    print('Selecting genes for model', model)
    if model in random_genes_dict:
        continue
    genes_shuffle = pd.Series()
    for group in range(1,7):
        gene_sub = genes_hall_3grps.Gene[genes_hall_3grps.Group==group]
        genes_shuffle = genes_shuffle.append(gene_sub.sample(n=len(gene_sub)//2))
    print('len shuffled genes', len(genes_shuffle))
    genes_dict, best_score = get_optimal_genes(data_3grps.loc[:, genes_shuffle],
                                               target_3grps,
                                               n_features=min_genes_3grps,
                                               big_step=5,
                                               step=5)
    best_genes = list(genes_dict[best_score[0]][1])
    random_genes_dict[model] = best_genes
    print(len(best_genes), 'genes selected')

# test with supervized model and select the best ones
model_scores_dict = {}
all_genes = []
good_models_dict = {}
for model in range(1, num_models + 1):
    print(model)
    cv_score = cross_val_score(SVC(kernel='linear'),
                               data_3grps.loc[:, random_genes_dict[model]],
                               target_3grps)
    print('Cross-val score:', cv_score)
    model_scores_dict[model] = cv_score
    all_genes = all_genes + random_genes_dict[model]
    if np.mean(cv_score) >= 0.95:
        good_models_dict[model] = random_genes_dict[model]

print('Number of genes selected:', len(set(all_genes)))
print('Number of random models selected', len(good_models_dict.keys()))


# apply 4-split supervised validation on samples
    # df_pred_dict where keys are the split # and values are dataframes with
    # class predictions
now = clean_name(str(datetime.now()))
df_pred_dict = {}
# split data and target into sets to loop over
idx_split = np.linspace(0, data_3grps.shape[0], 5).astype(int)
# loop over datasets
for split in np.arange(1, len(idx_split)):
    df_pred = pd.DataFrame(index=data_3grps.index)
    df_gene_set = pd.DataFrame(index=np.arange(1, 100))
    for k, v in good_models_dict.items():
        print(len(v))
        print('Set', k)
        gene_set = v
        data_psite = data_3grps.loc[:, gene_set]
        X_test = data_psite.iloc[
                idx_split[split - 1]: idx_split[split], :]
        X_train = data_psite.drop(X_test.index, axis=0)
        y_test = target_3grps.iloc[
                idx_split[split - 1]: idx_split[split]]
        y_train = target_3grps.drop(y_test.index, axis=0)
        y_test.name = 'True Class'
        
        clf = SVC(kernel='linear', probability=True)
        clf.fit(X_train, y_train)
        print('Testing score:', clf.score(X_test, y_test))
        y_pred = pd.Series(
                clf.predict(
                        X_test),
                        index=X_test.index)
        df_pred.loc[:, 'Gene Set #{}'.format(k)] = y_pred
        gene_v = [v[g] if g < len(v) else float('NaN')
                    for g in range(len(df_gene_set))]
        df_gene_set.loc[:, k] = gene_v
    percent_pred = df_pred.apply(
            lambda x: mode(x.values)[1][0] / df_pred.shape[1], axis=1)
    percent_pred.name = 'Percent Predicted'
    most_predicted = df_pred.mode(axis=1)
    most_predicted.columns = ['Most Predicted']
    df_pred = percent_pred.to_frame().join(df_pred)
    df_pred = most_predicted.join(df_pred)
    df_pred.index.name = 'sample_id'
    df_pred_dict[split] = df_pred

tot = []
av_len = 0
for col in df_gene_set.columns:
    av_len += len(df_gene_set.loc[:, col].dropna())
    tot += list(df_gene_set.loc[:, col].dropna())
tot = set(tot)
print(len(tot))
av_len = av_len / df_gene_set.shape[1]
print('Av length of genes', av_len)

writer = writer = pd.ExcelWriter(
        'Combined_validated_3grps_supervised_pred_{}_genes_{}.xlsx'.format(len(tot), now),
                            engine='xlsxwriter')
for k, v in df_pred_dict.items():
    print('Split', k)
    v.dropna().to_excel(writer,
               sheet_name='split {}'.format(k),
               index=True)
df_gene_set.to_excel(writer,
                      sheet_name='Genes Used',
                      index=False)
writer.save()

# get gene frequencey
gene_freq = []
for col in df_gene_set.columns:
    gene_freq += list(df_gene_set.loc[:, col].dropna())
gene_freq = pd.Series(gene_freq).value_counts()

# run models on ambiguous samples
# evaluate ambiguous cases
data_amb = data.loc[clin.DBU_4grp == 'Ambi', data_3grps.columns]
random_state = 0
df_pred_ambi = pd.DataFrame(index=data_amb.index)
df_proba = pd.DataFrame(np.zeros([data_amb.shape[0], 4]),
                        index=data_amb.index,
                        columns=['Total_Prob1', 'Total_Prob2', 'Total_Prob3','Total_Prob4'])
df_gene_set = pd.DataFrame(index=np.arange(1, 100))
for k, v in good_models_dict.items():
    random_state += 1
    print(len(v))
    print('Set', k)
    gene_set = v
    
    # split data
    X_train, X_test, y_train, y_test = train_test_split(
        data_3grps.loc[:, gene_set],
        target_3grps,
        random_state=random_state)
    
    # train classifier
    clf = SVC(kernel='linear', probability=True)
    clf.fit(X_train, y_train)
    print('Testing score:', clf.score(X_test, y_test))
    
    # get predictions
    y_pred = pd.Series(
            clf.predict(
                    data_amb.loc[:, gene_set]),
                    index=data_amb.index)
    
    df_pred_ambi.loc[:, 'Gene Set #{}'.format(k)] = y_pred
    
    # get probabilities
    y_proba = pd.DataFrame(clf.predict_proba(data_amb.loc[:, gene_set]),
                           index=data_amb.index,
                           columns=['Set{} Proba1'.format(k),
                                    'Set{} Proba2'.format(k),
                                    'Set{} Proba3'.format(k),
                                    'Set{} Proba4'.format(k)])
    df_proba = df_proba + y_proba.values
    
    # get genes used
    gene_v = [v[g] if g < len(v) else float('NaN')
                for g in range(len(df_gene_set))]
    df_gene_set.loc[:, k] = gene_v

percent_pred = df_pred_ambi.apply(
        lambda x: mode(x.values)[1][0] / df_pred_ambi.shape[1], axis=1)
percent_pred.name = 'Percent Predicted'
most_predicted = df_pred_ambi.mode(axis=1)
most_predicted.columns = ['Most Predicted']
df_pred_ambi = percent_pred.to_frame().join(df_pred_ambi)
df_pred_ambi = most_predicted.join(df_pred_ambi)
df_pred_ambi.index.name = 'sample_id'
df_proba = df_proba / len(good_models_dict.keys())
df_pred_ambi = df_pred_ambi.join(df_proba)
writer = writer = pd.ExcelWriter(
        'Combined_validated_ambiguous_supervised_pred_{}_genes_{}.xlsx'.format(len(tot), now),
                            engine='xlsxwriter')
df_pred_ambi.to_excel(writer,
                 sheet_name='Ambiguous Predictions',
                 index=True)
df_gene_set.to_excel(writer,
                      sheet_name='Genes Used',
                      index=False)
writer.save()

# run models on MCG samples
data_MCG = data.loc[clin.Study == 'MCG', data_3grps.columns]
random_state = 0
df_pred_MCG = pd.DataFrame(index=data_MCG.index)
df_proba = pd.DataFrame(np.zeros([data_MCG.shape[0], 4]),
                        index=data_MCG.index,
                        columns=['Total_Prob1', 'Total_Prob2', 'Total_Prob3','Total_Prob4'])
df_gene_set = pd.DataFrame(index=np.arange(1, 100))
for k, v in good_models_dict.items():
    random_state += 1
    print(len(v))
    print('Set', k)
    gene_set = v
    
    # split data
    X_train, X_test, y_train, y_test = train_test_split(
        data_3grps.loc[:, gene_set],
        target_3grps,
        random_state=random_state)
    
    # train classifier
    clf = SVC(kernel='linear', probability=True)
    clf.fit(X_train, y_train)
    print('Testing score:', clf.score(X_test, y_test))
    
    # get predictions
    y_pred = pd.Series(
            clf.predict(
                    data_MCG.loc[:, gene_set]),
                    index=data_MCG.index)
    
    df_pred_MCG.loc[:, 'Gene Set #{}'.format(k)] = y_pred
    
    # get probabilities
    y_proba = pd.DataFrame(clf.predict_proba(data_MCG.loc[:, gene_set]),
                           index=data_MCG.index,
                           columns=['Set{} Proba1'.format(k),
                                    'Set{} Proba2'.format(k),
                                    'Set{} Proba3'.format(k),
                                    'Set{} Proba4'.format(k)])
    df_proba = df_proba + y_proba.values
    
    # get genes used
    gene_v = [v[g] if g < len(v) else float('NaN')
                for g in range(len(df_gene_set))]
    df_gene_set.loc[:, k] = gene_v

percent_pred = df_pred_MCG.apply(
        lambda x: mode(x.values)[1][0] / df_pred_MCG.shape[1], axis=1)
percent_pred.name = 'Percent Predicted'
most_predicted = df_pred_MCG.mode(axis=1)
most_predicted.columns = ['Most Predicted']
df_pred_MCG = percent_pred.to_frame().join(df_pred_MCG)
df_pred_MCG = most_predicted.join(df_pred_MCG)
df_pred_MCG.index.name = 'sample_id'
df_proba = df_proba / len(good_models_dict.keys())
df_pred_MCG = df_pred_MCG.join(df_proba)
writer = writer = pd.ExcelWriter(
        'Combined_validated_MCG_supervised_pred_{}_genes_{}.xlsx'.format(len(tot), now),
                            engine='xlsxwriter')
df_pred_MCG.to_excel(writer,
                 sheet_name='MCG Predictions',
                 index=True)
df_gene_set.to_excel(writer,
                      sheet_name='Genes Used',
                      index=False)
writer.save()