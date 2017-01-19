#!/usr/bin/env python3

import pandas as pd

class ReadCountTable(pd.DataFrame):
    
    def __init__(self, df, gene_column_name = 'genename'):
        
        super().__init__(df)
        self.__gene_column_name = gene_column_name
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
            if columnname == self.__gene_column_name:
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
        
        
    def write_DESeq2_files(self, directory, rct_file_name = 'read_count_table.txt', coldata_file_name = 'col_data.txt'):
        
        ''' Writes out necessary files to perform a DESeq2 analysis on the read count table.
        Two files are generated:
        1) the read count table
        2) a table containing information on the columns and how they're grouped.
        Files are written to the given directory. '''
        
        # write count table
        with open('{}/{}'.format(directory, rct_file_name), 'w') as f:
            f.write('# Read counts\n')
            self.to_csv(f, sep='\t', index=False)
        
        # write column data table
        with open('{}/{}'.format(directory, coldata_file_name), 'w') as f:
            f.write('Ind\tDay\tRepeat\tRun\n')
            for columnname in self.columns:
                if columnname != 'genename':
                    ind = columnname.split('/')[-1].split('_')[0]
                    day = columnname.split('/')[-1].split('_')[1].replace('EXP', '')
                    rep = columnname.split('/')[-1].split('_')[2].replace('EXP', '')
                    run = columnname.split('/')[-3].split('_')[0]
                    f.write('{}\t{}\t{}\t{}\n'.format(ind, day, rep, run))
        
        
def concat_read_count_tables(listoftables):
    
    ''' Given multiple read count table instances, combine them into a single instance.
    List of genes in the final table contains all the genes present in different read count tables.
    If a read count table didn't have any data on a gene present in another table,
    it's samples are each assigned a value of 0 for that gene. '''
    
    # concat all rcts
    total_rct = listoftables[0].set_index('genename')
    for i in range(len(listoftables)-1):
        rct_to_add = listoftables[i+1].set_index('genename')
        total_rct = pd.concat([total_rct, rct_to_add], axis=1)
        
    # replace missing values with 0     
    total_rct = total_rct.fillna(0)
    
    # set genename as a column instead of index
    total_rct['genename'] = total_rct.index

    return ReadCountTable(total_rct)


if __name__ == '__main__':
    
    from glob import glob
    
    indir = "/Users/nicolasdeneuter/Bestanden/PhD/Projects/GOA/RNAseq/readcounts"
    outdir = "/Users/nicolasdeneuter/Bestanden/PhD/Projects/GOA/RNAseq/readcounts"
    
    rct_list = []
    for filepath in glob(indir+'/*'):
        if 'readcount' in filepath.split('/')[-1]:
            rct = ReadCountTable(pd.read_csv(filepath, sep = '\t'))
            rct.combine_lane_counts()
            rct_list.append(rct)
    
    total_rct = concat_read_count_tables(rct_list)
    total_rct.remove_sample('H1_EXP0_1')
    
    # write DESeq colData file
    total_rct.write_DESeq2_files(outdir)