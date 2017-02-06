#!/usr/bin/env python3

import os
import sys
import subprocess
import pandas as pd
import readcountclass
from glob import glob

# read pipeline argument
pipeline = sys.argv[1]
    
# determine pipeline to be used based on argument
if pipeline == '-R' or pipeline == '-R0':
    print('Running R based pipeline')
elif pipeline == '-ML':
    print('Running ML based pipeline')
else:
    raise ValueError('Expects either -R, -R0 or -ML as argument to determine which pipeline to run.')

# read data from working directory
maindir = os.getcwd()
print('Looking for files in {}'.format(maindir))
#maindir = "/Users/nicolasdeneuter/Bestanden/PhD/Projects/GOA/RNAseq/readcounts"

# read each count table in maindir and save it as count table instance
rct_list = []
print('Files found:')
for filepath in glob(maindir+'/*'):
    if 'readcount' in filepath.split('/')[-1]:
        print(filepath)
        rct = readcountclass.ReadCountTable(pd.read_csv(filepath, sep = '\t'))
        # combine counts for same sample over different lanes
        rct.combine_lane_counts()
        rct_list.append(rct)

if len(rct_list) == 0:
    raise IOError('No read count files were found.')

# combine all count tables
print('Combining read count tables together')
total_rct = readcountclass.concat_read_count_tables(rct_list)

print('Removing samples with unknown responder status')
response_dict = readcountclass.make_responders_dict()
# remove samples for which responder status is unknown
for sample in total_rct.columns:
    if sample.split('/')[-1].split('_')[0] not in response_dict.keys() and sample != total_rct.gene_column_name:
        del total_rct[sample]

print('Combining data on H6_EXP3_1 (sequenced twice with bad quality)')
# H6_EXP3_1 was sequenced twice due to bad quality, take mean of the two runs since they're both of lower quality
H6_EXP3_1_samples = ['/home/shared_data_immuno/Run_Data/170111_NB501809_0047_AH2HH2BGX2/tmp_files/H6_EXP3_1_S29',\
                '/home/shared_data_immuno/Run_Data/161125_NB501809_0023_AHLC7CBGXY/tmp_files/H6_EXP3_1_S15']
total_rct['/home/shared_data_immuno/Run_Data/161125+170111_NB501809_0047_AH2HH2BGX2/tmp_files/H6_EXP3_1'] = \
    total_rct[H6_EXP3_1_samples].sum(axis =1).map(lambda x: int(x/2))
total_rct = readcountclass.ReadCountTable(total_rct.drop(H6_EXP3_1_samples, axis = 1))


if pipeline == '-R' or pipeline == '-R0':
    
    # write data to files in a subdir
    if pipeline == '-R':
        out = '{}/R_pipeline_nonresp_day'.format(maindir)
    elif pipeline == '-R0':
        out = '{}/R_pipeline_nonresp_0'.format(maindir)
        
    print('Data preprocessing finished.\Placing output in directory: {}'.format(out))
    if os.path.isdir(out) != True:
        os.mkdir(out)
    total_rct.write_DESeq2_files(out, response_dict, pipeline = pipeline)
    
    # perform R analysis
    os.chdir(out)
    
    ## Need to pass arguments to R scripts to use different model based on R pipeline ##
    bashcommand = 'Rscript {}/deseq2_R_analysis.R'.format(maindir)
    process = subprocess.run(bashcommand.split())
    
    
if pipeline == '-ML':
    
    # create dict mapping volunteer ids to their samples
    volunteer_map = {}
    for columnname in total_rct.columns:
        if columnname != total_rct.gene_column_name:
            volunteer_id = columnname.split('/')[-1].split('_')[0]
            volunteer_map.setdefault(volunteer_id, []).append(columnname)
            
    # make directory to place all output in
    out = '{}/ML_pipeline'.format(maindir)
    print('Data preprocessing finished.\Placing output in directory: {}'.format(out))
    if os.path.isdir(out) != True:
        os.mkdir(out)
    os.chdir(out)
    
    # write away each volunteer's data to a separate subdir
    for volunteer_id, samples in volunteer_map.items():
        subout = '{}/{}'.format(out, volunteer_id)
        if os.path.isdir(subout) != True:
            os.mkdir(subout)
        vol_rct = readcountclass.ReadCountTable(total_rct[[total_rct.gene_column_name]+samples])
        vol_rct.write_DESeq2_files(subout, response_dict, pipeline = pipeline)
        
        # perform R analysis on this part of the data
        os.chdir(subout)
        bashcommand = 'Rscript {}/deseq2_ML_analysis.R'.format(maindir)
        process = subprocess.run(bashcommand.split())        
        os.chdir(out)
        
    
        
print('Finished Python script')

