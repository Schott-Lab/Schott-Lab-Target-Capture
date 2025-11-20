# Masked Coverage
from re import T
import zipfile
import pandas as pd
import numpy as np
import os
import sys
import time
import pandas as pd
#!pip install xlwt
#!pip install xlrd
#!pip install xlutils
import xlwt
from xlwt.Workbook import *
from pandas import ExcelWriter
#!pip install xlsxwriter

#Get all .tsv file names from command line arguments
names_foo = sys.argv[1:]
df = pd.DataFrame(names_foo)
df[0] = df[0].str.replace('.*/', '', regex=True)
#Initialize empty lists to store dataframes and overall coverage
exp = []
ovr = []

#Read all .tsv files into a list of dataframes stored in exp
for i in range(len(names_foo)):
  h = pd.read_csv(names_foo[i], sep = '\t')
  exp.append(h)

#Create dataframe to store gene names
### Changed Code here ###
#gene_names = exp[0].iloc[:-1, 0] #Source of error, gets gene names from first .tsv file

gene_names = None

#Extract gene names from each dataframe and merge into a single dataframe to ensure all gene names are included
for species in exp:
    species_gene_name = species.iloc[:-1, 0]
    species_gene_name = pd.DataFrame(species_gene_name)
    species_gene_name['Seq_ID'] = species_gene_name['Seq_ID'].str.replace('.*_', '', regex=True)
    species_gene_name = species_gene_name.rename(columns={"Seq_ID": "Gene"})
    if gene_names is None:
        gene_names = species_gene_name  # Initialize gene with the first species' gene names
    gene_names = pd.merge(gene_names,species_gene_name, how="outer")
### Changed Code here ###

#Create dataframe to store species names
species = pd.DataFrame(names_foo)
species.columns = ['Species']
species['Species'] = species['Species'].str.replace('.*coverage_', '', regex=True)
species['Species'] = species['Species'].str.replace('_trimmed.tsv', '', regex=True)

#Calculate total and average values for each .tsv file and add to dataframe and storing it in l3sum and l3ave

### Can be improved ###
l1 = []
l2 = list(np.zeros(len(exp)))
l3 = list(np.zeros(len(exp)))
l3sum = np.zeros(len(exp))
l3ave = np.zeros(len(exp))
#Calculate Length sums and averages of each dataframe
for i in range(len(exp)):
  l1.append(exp[i]['Length'][0:len(exp[i]['Length'])-1])
  l2[i] = l1[i]
  l3[i] = [float(x) for x in l2[i]]
  l3sum[i] = np.sum(l3[i])
  l3ave[i] = np.mean(l3[i])
### Can be improved ###


#Initialize arrays to store sums and averages of each column in each dataframe
nsum = np.zeros(len(exp))
covsum = np.zeros(len(exp))
npsum = np.zeros(len(exp))
cpsum = np.zeros(len(exp))

nave = np.zeros(len(exp))
covave = np.zeros(len(exp))
npave = np.zeros(len(exp))
cpave = np.zeros(len(exp))

#Calculate sums and averages for each column in each dataframe
for t in range(len(exp)):

  nsum[t] = np.sum(exp[t]['N_count'][0:len(exp[t]['N_count'])-1])
  nave[t] = np.mean(exp[t]['N_count'][0:len(exp[t]['N_count'])-1])

  covsum[t] = np.sum(exp[t]['cov_count'][0:len(exp[t]['cov_count'])-1])
  covave[t] = np.mean(exp[t]['cov_count'][0:len(exp[t]['cov_count'])-1])

  npsum[t] = np.sum(exp[t]['N_perc'][0:len(exp[t]['N_perc'])-1])
  npave[t] = np.mean(exp[t]['N_perc'][0:len(exp[t]['N_perc'])-1])

  cpsum[t] = np.sum(exp[t]['Cov_perc'][0:len(exp[t]['Cov_perc'])-1])
  cpave[t] = np.mean(exp[t]['Cov_perc'][0:len(exp[t]['Cov_perc'])-1])


#Add total and average rows to each dataframe and store overall average coverage
### Change Code here ###
gene = gene_names

for m in range(len(exp)):
  #Add Total and Average rows to each dataframe
  new_row = ['Total',l3sum[m],nsum[m],covsum[m],npsum[m],cpsum[m]]
  new_row2 = ['Average',l3ave[m],nave[m],covave[m],npave[m],cpave[m]]
  exp[m].loc["Total"] = new_row
  exp[m].loc["Average"] = new_row2

  #Store overall average coverage
  ovr.append(cpave[m])
  #Add coverage percentage values to gene dataframe
  perc_val = list(exp[m]['Cov_perc'])
  del perc_val[-3:] #Remove last three entries which are Average Coverage, Total, and Average
  speciesName = species.loc[m]['Species']

  gene_names_for_species = list(exp[m].iloc[:-1, 0].str.replace('.*_', '', regex=True))
  del gene_names_for_species[-2:]

  #Indexes perc_val of each species name to the correct gene in the gene dataframe
  #Matches the perc_val of each species to it's own corresponding gene name
  #Merges the two dataframes on the 'gene dataframe to align perc_val to the genes in the gene dataframe
  #and get the indexed column to add to the gene dataframe
  perc_val_gene_match = pd.DataFrame({'Gene': gene_names_for_species, 'Perc_Val': perc_val})
  perc_val_gene_match = perc_val_gene_match.set_index('Gene')
  perc_val_indexed = pd.merge(gene, perc_val_gene_match, on="Gene", how="left")['Perc_Val']
  
  gene[speciesName] = perc_val_indexed
  
### Changed Code here ###

gene = gene.set_index(['Gene'])
gene['Gene Mean'] = gene.mean(axis=1)
gene.loc['Species Mean'] = gene.mean()
gene = gene.fillna(0)

gene.to_csv("Coverage_Gene_Summary.csv", index=True)

Ovr_df = df.join(pd.DataFrame({'Average Coverage': ovr}))
ovrave = ['Overall Average', Ovr_df.loc[:, 'Average Coverage'].mean()]
Ovr_df.loc["Overall Average"] = ovrave

sheetname = []
for m in range(len(exp)):
  sheetname.append(df[0][m])

c = list(np.zeros(len(exp)))
for s in range(len(exp)):
  c[s] = df[0][s]

for i in range(len(exp)):
    exp[i].to_csv("%s.csv" % c[i], index=False)

Ovr_df.to_csv("Coverage_Species_Summary.csv", index=False)