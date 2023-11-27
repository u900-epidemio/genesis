#!/usr/bin/env python

### Import packages ###
import numpy as np
import pandas as pd
import sys


### Import file with the type of variants ###
df_mpileup = pd.read_csv(sys.argv[len(sys.argv)-2], sep=",", header = None)

gene_position = pd.read_csv(sys.argv[len(sys.argv)-1], sep=",")


### Define sample name ###
n = sys.argv[1].find('.')
name = sys.argv[1][0:n]


### Change colums name in df_mpileup ###
df_mpileup.columns = ['chromosome', 'position', 'base', 'reads_number']


### Assign gene ID for each position ###
df_mpileup['gene'] = ''
counter1 = np.zeros(len(gene_position)) # Number of bases in the gene
counter2 = np.zeros(len(gene_position)) # Number of bases with (reads_number >= 30) in the gene

for i in range(len(df_mpileup)): 
    for j in range(len(gene_position)):
        if ((gene_position['chr'][j] == df_mpileup['chromosome'][i]) and (gene_position['end'][j] >= df_mpileup['position'][i] >= gene_position['start'][j])):
            df_mpileup['gene'][i] = gene_position['gene_ID'][j]
            counter1[j] += 1
            if df_mpileup['reads_number'][i] >= 30 :
                counter2[j] += 1
    if (df_mpileup['gene'][i] == ''):
        df_mpileup['gene'][i] = 'not in the data base'


### Calculate the coverage ###
gene_position['coverage'] = ''
for j in range(len(gene_position)):
    gene_position['coverage'][j] =  (100 * (counter2[j] / counter1[j])).round(decimals = 2)


### Export in csv ###
gene_position.to_csv(f"{name}.coverage.mpileup.python.csv", index = False)
