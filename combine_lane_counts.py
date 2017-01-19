#!/usr/bin/env python

import glob
import sys
import natsort
import pandas as pd

# directories to find data in and output results to
indir = "/Users/nicolasdeneuter/Bestanden/PhD/Projects/GOA/RNAseq/readcounts"
outdir = "/Users/nicolasdeneuter/Bestanden/PhD/Projects/GOA/RNAseq/readcounts"

# number of patients to include
target_pats = 1

#######

count_dfs = []

# for each readcount file
for filepath in glob.glob(indir+'/*'):
    
    # only use files containing read count data
    if 'readcounts' not in filepath.split('/')[-1]:
        continue
        
    print('Reading file {}'.format(filepath))
    
    # read in data
    data = pd.read_csv(filepath, sep = '\t')
    
    # find columnnames that contain data for same patient & time point but have been sequenced on different lanes
    lanegroups = {}
    for columnname in data.columns:
        
        # skip irrelevant columns for grouping
        if columnname == 'genename':
            continue
        elif 'Unnamed' in columnname:
            print(data[columnname])
            del data[columnname]
            continue
        # make common name for lanegroup
        sample = '_'.join(columnname.split('_')[:-3])
        # make dict to keep track of which columns fall in one group
        lanegroups.setdefault(sample, []).append(columnname)
        
    # create new dataframe with gene names
    df = pd.DataFrame({'genename':data['genename']})
    
    # sum read counts for lanes from one group together and add them as a column to the new df
    for columngroup, columnnames in lanegroups.items():
        df[columngroup] = data[columnnames].sum(axis=1)
    
    print('Genes found: {}'.format(df.shape[0]-1))
    print('Combined {} lane-specific columns into {}'.format(len(data.columns)-1,len(df.columns)-1))
        
    # keep df for this file in a list to concatenate all of them later on
    count_dfs.append(df)

print('Combining data from all files.')

# concat all dfs together, using genenames as index
total_df = count_dfs[0].set_index('genename')
for i in range(len(count_dfs)-1):
    df_to_add = count_dfs[i+1].set_index('genename')
    total_df = pd.concat([total_df, df_to_add], axis=1)

# replace missing values with 0     
total_df = total_df.fillna(0)

# set genename as a column instead of index
total_df['genename'] = total_df.index

# rearrange columns to put 'genename' as the first column
resorted_cols = ['genename']+natsort.natsorted(list(total_df.columns[:-1]), key=lambda x:x.split('/')[-1])
resorted_cols.remove('/home/shared_data_immuno/Run_Data/161125_NB501809_0023_AHLC7CBGXY/tmp_files/H6_EXP3_1_S15')

# depending on the number of patients to be included, make a shorter list with only the relevant sample names
reduced_cols = ['genename']
num_of_pats = 0
i = 1
while num_of_pats < target_pats and num_of_pats <= len(resorted_cols):
    reduced_cols.append(resorted_cols[i])
    num_of_pats = (len(reduced_cols)-1)/6
    i += 1
        
print('Total number of patients included: {}'.format((len(reduced_cols)-1)/6))

# rearrange dataframe to only contain requested data
total_df = total_df[reduced_cols]

print('Total number of genes found: {}'.format(total_df.shape[0]-1))

# write count data to file
print('Writing data to {}/combined_lane_counts.txt'.format(outdir))
with open('{}/combined_lane_counts.txt'.format(outdir), 'w') as f:
    f.write('# Read counts combined over different lanes\n')
    total_df.to_csv(f, sep='\t', index=False)
    
# write DESeq colData file
print('Writing colData for DESeq2 analysis in R to {}/colData.txt'.format(outdir))
with open('{}/colData.txt'.format(outdir), 'w') as f:
    f.write('Ind\tDay\tRepeat\tRun\n')
    for columnname in reduced_cols:
        if columnname != 'genename':
            ind = columnname.split('/')[-1].split('_')[0]
            day = columnname.split('/')[-1].split('_')[1].replace('EXP', '')
            rep = columnname.split('/')[-1].split('_')[2].replace('EXP', '')
            run = columnname.split('/')[-3].split('_')[0]
            f.write('{}\t{}\t{}\t{}\n'.format(ind, day, rep, run))
        
