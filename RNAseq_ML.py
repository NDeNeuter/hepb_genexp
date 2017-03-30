#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
from readcountclass import DiffGeneTable, make_responders_dict
from sklearn.decomposition import PCA

#main_data_dir = sys.argv[1]
main_data_dir = "/Users/nicolasdeneuter/Bestanden/PhD/Projects/GOA/RNAseq/readcounts/ML_pipeline"

resp_dict = make_responders_dict(threshold = 50)

subdirmap = {}

# for each subdir in the main directory, find result files (each contrast has a different result file)
for x in os.walk(main_data_dir):
    subdir = x[0]
    if subdir == main_data_dir:
        continue
    data_files = glob(subdir+'/results*')
    # for each contrast's result:
    # make a dataframe containing a single row (ie one volunteer) containing log fold changes for each column/gene
    for data_file in data_files:
        # determine characteristics of file
        volunteer = data_file.split('/')[-2]
        contrast = data_file.split('/')[-1].replace('.txt','').replace('results_','')
        resp = resp_dict[volunteer]
        # # skip day 3 vs 7 contrast because it's the difference of the other contrasts
        # if contrast == '3_7':
        #     continue
        # read in and format data
        diff_gene_table = pd.read_csv(data_file).set_index('Unnamed: 0')
        # # only use those genes that are significantly differentially expressed
        # diff_gene_table = diff_gene_table[diff_gene_table['padj'] <= 0.05]
        gene_folds = diff_gene_table.transpose().loc[['log2FoldChange']]
        gene_folds = gene_folds.rename(index={'log2FoldChange': '{}'.format(volunteer)}, 
                                       columns={name: '{}_{}'.format(name, contrast) for name in gene_folds.columns})
        gene_folds.columns.name = 'volunteer_resp'
        subdirmap.setdefault(volunteer, []).append(gene_folds)
        
temp_fold_frames = []
# combine dataframe on same volunteer but different contrasts into one dataframe
# result for each volunteer is a new dataframe with one row and columns containing gene+contrast info
for volunteer, gene_fold_list in subdirmap.items():
    temp_fold = pd.concat(gene_fold_list, axis = 1)
    temp_fold['response'] = resp_dict[volunteer]
    temp_fold_frames.append(temp_fold)
    
gene_fold = temp_fold_frames[0]
for i in range(1, len(temp_fold_frames)):
    gene_fold_to_add = temp_fold_frames[i]
    gene_fold = gene_fold.append(gene_fold_to_add)
    
# replace missing values with 0
total_gene_fold = gene_fold.fillna(0)

X = total_gene_fold[[x for x in total_gene_fold.columns if x != 'response']]
y = total_gene_fold['response']

pca = PCA(n_components = 2)
X_r = pd.DataFrame(pca.fit(X).transform(X), index = X.index)

plt.figure()
colors = ['turquoise', 'darkorange']
target_names = ['Non-resp', 'Resp']
label = [float(x.replace('H','')) for x in X_r.index]
for color, target_name in zip(colors, target_names):
    plt.scatter(X_r[y == target_name][0], X_r[y == target_name][1], color=color, alpha=.8, lw=2, label=target_name)
plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.xlabel(pca.explained_variance_ratio_[0])
plt.ylabel(pca.explained_variance_ratio_[1])
plt.title('PCA of gene-expression data')
for i, txt in enumerate(label):
    plt.annotate('H'+str(int(txt)), (X_r[0][i],X_r[1][i]))
plt.show()