#!/usr/bin/env python

import glob
import pandas as pd

# directories to find data in and output results to
indir = "/Users/nicolasdeneuter/Bestanden/PhD/Projects/GOA/RNAseq/readcounts/"
outdir = "/Users/nicolasdeneuter/Bestanden/PhD/Projects/GOA/RNAseq/readcounts/"

count_dfs = []

# for each readcount file
for filepath in glob.glob(indir+'/*'):
    
    # only use files containing read count data
    if 'readcounts' not in filepath.split('/')[-1]:
        continue
        
    print('Processing file {}'.format(filepath))
    
    # read in data
    data = pd.read_csv(filepath, sep = '\t')
    
    # find columnnames that contain data for same patient & time point but have been sequenced on different lanes
    lanegroups = {}
    for columnname in data.columns:
        
        # skip irrelevant columns for grouping
        if columnname == 'genename' or 'Unnamed' in columnname:
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
        
    # keep dfs in a list to concat them later on
    count_dfs.append(df)

# concat all dfs together, using genenames as index
total_df = count_dfs[0].set_index('genename')
for i in range(len(count_dfs)-1):
    df_to_add = count_dfs[i+1].set_index('genename')
    total_df = pd.concat([total_df, df_to_add], axis=1)

# add genename as a column
total_df['genename'] = total_df.index

# rearrange columns to put 'genename' as the first column
resorted_cols = ['genename']+sorted(total_df.columns[:-1])
total_df = total_df[resorted_cols]

# write to file
with open('{}/combined_lane_counts.txt'.format(outdir), 'w') as f:
    f.write('# Read counts combined over different lanes\n')
    total_df.to_csv(f, sep='\t', index=False)
    