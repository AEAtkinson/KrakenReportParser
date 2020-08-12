#!/usr/bin/env python

# Gather kraken2 reports from a directory and assemble them into 1 large table
# No filters applied
# Rows are observations from the kraken2 report
# Columns are individuals
# 0 is used where a taxon was not observed in a given sample

# USAGE: collect_kraken2_reports.py --dir [directory to look for reports in] --out [results file name] 


import pandas as pd
import glob, os

# build arg parser


# glob in files
report_list = glob.glob('*.txt')

samples = []
# loop through kraken2 reports
for report in report_list:
    
    # get sample ID from filename
    sample_name = report.split('.')[0]
    print(sample_name)
    # read in as pandas dataframe
    result = pd.read_csv(report, sep="\t", header=None)
    result.columns = ["percent_total", "reads", "taxon_reads", "taxon", "NCBI_taxon_ID","name"]
    
    # drop undesired columns and strip white space
    result = result[['name','reads']]
    result.name = result.name.str.strip()
    
    ##sample.columns = ['name', sample_name]
    # strip white space from name columns
    ##sample.name = sample.name.str.strip()
    result = result[result['name'].isin(['Human papillomavirus type 11', 'Human papillomavirus type 6b', 'Human papillomavirus type 31', 'Human papillomavirus type 33', 'Human papillomavirus type 58', 'Human papillomavirus type 52', 'Human papillomavirus type 16', 'Human papillomavirus type 39', 'Human papillomavirus type 45', 'Human papillomavirus type 59', 'Human papillomavirus type 18', 'Human papillomavirus type 51', 'Human papillomavirus type 56', 'Human T-cell leukemia virus type I', 'Human herpesvirus 4 type 2', 'Human herpesvirus 4 type 1', 'Human herpesvirus 5 (strain 1042)', 'Human polyomavirus 1', 'HBV genotype A'])]
    #sort columns
    result.sort_values(by=['name'])	
    #sample.columns = ['name', sample_name]
    #print(result)
    #rotate table
    df = result.transpose()
    #reset columns
    
    
    df.columns = df.iloc[0]
    #df.rename(columns=df.iloc[0]).drop(df.index[0])
    df['Sample'] = sample_name
    #df = df.reset_index()
    df2 = df[1:]
    #order columns
    df2 = df2[['Sample', 'Human papillomavirus type 6b', 'Human papillomavirus type 11', 'Human papillomavirus type 16', 'Human papillomavirus type 18', 'Human papillomavirus type 31', 'Human papillomavirus type 33', 'Human papillomavirus type 39', 'Human papillomavirus type 45', 'Human papillomavirus type 51', 'Human papillomavirus type 52', 'Human papillomavirus type 56', 'Human papillomavirus type 58', 'Human papillomavirus type 59', 'HBV genotype A', 'Human polyomavirus 1', 'Human herpesvirus 5 (strain 1042)', 'Human herpesvirus 4 type 2', 'Human herpesvirus 4 type 1', 'Human T-cell leukemia virus type I']]
    df2.rename(columns={'Human papillomavirus type 6b': 'HPV6', 'Human papillomavirus type 11': 'HPV11', 'Human papillomavirus type 16': 'HPV16', 'Human papillomavirus type 18': 'HPV18', 'Human papillomavirus type 31': 'HPV31', 'Human papillomavirus type 33': 'HPV33', 'Human papillomavirus type 39': 'HPV39', 'Human papillomavirus type 45': 'HPV45', 'Human papillomavirus type 51': 'HPV51', 'Human papillomavirus type 52': 'HPV52', 'Human papillomavirus type 56': 'HPV56', 'Human papillomavirus type 58': 'HPV58', 'Human papillomavirus type 59': 'HPV59', 'HBV genotype A': 'HBV A', 'Human polyomavirus 1': 'BKV', 'Human herpesvirus 5 (strain 1042)': 'CMV', 'Human herpesvirus 4 type 2': 'EBV-2', 'Human herpesvirus 4 type 1': 'EBV-1', 'Human T-cell leukemia virus type I': 'HTLV1'}, inplace=True)
    #print(df2)
    #save to file
    csv_file = os.path.splitext(report)[0]+".temp.tsv"
    #sample.to_csv(csv_file, sep='\t')
    df2.to_csv(csv_file, sep='\t', index=False)

filenames = glob.glob('*.temp.tsv')
dfs = []
for filename in filenames:
    dfs.append(pd.read_table(filename))
results = pd.concat(dfs, sort=False)
#print(results)
results.to_csv('HNCViralMasterV2.tsv', sep='\t', index=False)
