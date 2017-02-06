#!/usr/bin/env python3

import natsort
import subprocess
import pandas as pd
from numpy import mean

class ReadCountTable(pd.DataFrame):

    def __init__(self, df, gene_column_name = 'genename'):
        
        super().__init__(df)
        self.gene_column_name = gene_column_name
        for columnname in self.columns:
            if 'Unnamed' in columnname:
                del self[columnname]
      
      
    def combine_counts(self, columndict):
        
        ''' Combines columns in the read count table by summing them together into a new column.
        Requires a dictionary columndict with:
        - as keys: names of the new columns
        - as value for a key: list of columnnames to be combined together.
        Combined columns are removed from the read count table. '''
        
        for columngroup, columnnames in columndict.items():
            self[columngroup] = self[columnnames].sum(axis=1)
            for column in columnnames:
                del self[column]
    
    
    def combine_lane_counts(self):
    
        ''' Combines read count data from data columns for the same sample but from different lanes together.
        This method is an extension of the combine_counts method and assumes lanes are encoded
        in samples names following a certaing format. '''
        
        # find columnnames that contain data for same patient & time point but have been sequenced on different lanes
        lanegroups = {}
        for columnname in self.columns:
            
            # skip irrelevant columns for grouping
            if columnname == self.gene_column_name:
                continue
            # make common name for lanegroup
            sample = '_'.join(columnname.split('_')[:-3])
            # make dict to keep track of which columns fall in one group
            lanegroups.setdefault(sample, []).append(columnname)

        self.combine_counts(lanegroups)
    
    
    def remove_sample(self, query):
        
        ''' Removes any columns in the read count data that contain the query. '''
        
        for col in self.columns:
            if query in col:
                del self[col]
        
        
    def write_DESeq2_files(self, directory, response_dict, rct_file_name = 'read_count_table.txt', coldata_file_name = 'col_data.txt', pipeline = 'ML'):
        
        ''' Writes out necessary files to perform a DESeq2 analysis on the read count table.
        Two files are generated:
        1) the read count table - rct_file_name
        2) a table containing information on the columns and how they're grouped - coldata_file_name
        Also uses a response_dict to write an extra factor (responder vs non-responder).
        Files are written to the given directory.
        Arguments:
        1) directory: dir in which output files will be placed
        2) rct_file_name: name of read count table file that will be created
        3) coldata_file_name: name of coldata table file (contains factors) that will be created
        4) processor: function that processes sample names to extract factors from it; is specific to analysis/dataset. '''
                
        # put gene name column first! DESeq2 analysis cares about order of columns
        resorted_cols = [self.gene_column_name]+natsort.natsorted(list([x for x in self.columns if x != self.gene_column_name]), key=lambda x:x.split('/')[-1])
        
        self = ReadCountTable(self[resorted_cols])
        
        # write count table
        print('Read count table: {}/{}'.format(directory, rct_file_name))
        with open('{}/{}'.format(directory, rct_file_name), 'w') as f:
            f.write('# Read counts\n')
            self.to_csv(f, sep='\t', index=False)
        
        # write column data table
        print('Column data table: {}/{}'.format(directory, coldata_file_name))
        with open('{}/{}'.format(directory, coldata_file_name), 'w') as f:
            
            f.write('Ind\tDay\tRepeat\tRun\tResp\n')
            for columnname in self.columns:
                if columnname != self.gene_column_name:
                    ind, day, rep, run, resp = process_samplename(columnname, response_dict, pipeline)
                    f.write('{}\t{}\t{}\t{}\t{}\n'.format(ind, day, rep, run, resp))


class DiffGeneTable(pd.DataFrame):

    def __init__(self, df):
        
        super().__init__(df)
        for columnname in self.columns:
            if 'Unnamed' in columnname:
                del self[columnname]
                
    def diff_genes(self, p_adj_threshold = 0.05):
        
        ''' Returns a list of differentially expressed genes with an adjsted p-value under a given threshold. '''
        
        return list(self.sort_values(by='padj')[self['padj'] <= p_adj_threshold].index)
        
    
    def fraction_diff_genes(self, p_adj_threshold = 0.05):
        
        ''' Returns % of genes that are differentially expressed among all genes in the table. '''
        
        all_genes = len(self.index)
        diff_genes = len(self[self['padj'] <= p_adj_threshold].index)
        
        return diff_genes/all_genes
    
        
def concat_read_count_tables(listoftables):
    
    ''' Given multiple read count table instances, combine them into a single instance.
    List of genes in the final table contains all the genes present in different read count tables.
    If a read count table didn't have any data on a gene present in another table,
    it's samples are each assigned a value of 0 for that gene. '''
    
    # concat all rcts
    # always set column containing gene names as index
    total_rct = listoftables[0].set_index(listoftables[0].gene_column_name)
    for i in range(len(listoftables)-1):
        rct_to_add = listoftables[i+1].set_index(listoftables[i+1].gene_column_name)
        total_rct = pd.concat([total_rct, rct_to_add], axis=1)
    
    # replace missing values with 0
    total_rct = ReadCountTable(total_rct.fillna(0))
    
    # set genename as a column instead of index
    total_rct[total_rct.gene_column_name] = total_rct.index

    return ReadCountTable(total_rct)


def process_samplename(columnname, response_dict, pipeline = 'ML'):
    
    ''' Processes a sample's name into the different factors associated with the sample.
    To be used during the writing of coldata (ie: factor table) with the
    write_DESeq2_files method of the ReadCountTable class. '''

    ind = columnname.split('/')[-1].split('_')[0]
    resp = response_dict[ind]
    repeat = columnname.split('/')[-1].split('_')[2]
    run = columnname.split('/')[-3].split('_')[0]
    
    if pipeline == '-R0' and resp == 'Non-resp':
        day = 0
    else:
        day = columnname.split('/')[-1].split('_')[1].replace('EXP', '')
        
    return ind, day, repat, run, resp


def make_responders_dict(datafile = "/Users/nicolasdeneuter/Dropbox/GOA/HepB run/Hep B run incl test _ Samples overview.xlsx"):
    
    ''' Process data on volunteers to determine if they're responders or not.
    Return a dict with volunteers as keys and their response as key. '''
    
    # read in data
    all_data = pd.read_excel(datafile, skiprows = [2, 3, 4, 5], header = [0, 1], sheetname = 'ELISA', index_col = 0)
    # select subset of data on elisa results for day 60
    elisa_data = all_data['Anti-HBS (IU/L)'][['ELISA_60 (Batch)', 'ELISA_60 (Repeat)']]
    # drop patients with missing values
    elisa_data = elisa_data.dropna()
    # change everything to numbers (< 2 => 0 and > 1000 => 1000)
    elisa_data = elisa_data.applymap(lambda x: x if x != '< 2' else 0)
    elisa_data = elisa_data.applymap(lambda x: x if x != '>1000' else 1000)
    # set each patient to responder (True) or non-responder (False) and turn series to dict after calculating mean response
    elisa_data = dict(elisa_data.apply(mean, axis = 1).map(lambda x: 'Non-resp' if x <= 10 else 'Resp'))
    
    return elisa_data


if __name__ == '__main__':
    
    print(make_responders_dict())