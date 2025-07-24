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

names_foo = sys.argv[1:]
df = pd.DataFrame(names_foo)
df[0] = df[0].str.replace('.*/', '', regex=True)
exp = []
ovr = []

for i in range(len(names_foo)):
  h = pd.read_csv(names_foo[i], sep = '\t')
  exp.append(h)

gene = exp[0].iloc[:-1, 0]
gene = pd.DataFrame(gene)
gene['Seq_ID'] = gene['Seq_ID'].str.replace('.*_', '', regex=True)
gene = gene.rename(columns={"Seq_ID": "Gene"})

species = pd.DataFrame(names_foo)
species.columns = ['Species']
species['Species'] = species['Species'].str.replace('.*coverage_', '', regex=True)
species['Species'] = species['Species'].str.replace('_trimmed.tsv', '', regex=True)


l1 = []
l2 = list(np.zeros(len(exp)))
l3 = list(np.zeros(len(exp)))
l3sum = np.zeros(len(exp))
l3ave = np.zeros(len(exp))
for i in range(len(exp)):
  l1.append(exp[i]['Length'][0:len(exp[i]['Length'])-1])
  l2[i] = l1[i]
  l3[i] = [float(x) for x in l2[i]]
  l3sum[i] = np.sum(l3[i])
  l3ave[i] = np.mean(l3[i])

nsum = np.zeros(len(exp))
covsum = np.zeros(len(exp))
npsum = np.zeros(len(exp))
cpsum = np.zeros(len(exp))

nave = np.zeros(len(exp))
covave = np.zeros(len(exp))
npave = np.zeros(len(exp))
cpave = np.zeros(len(exp))

for t in range(len(exp)):

  nsum[t] = np.sum(exp[t]['N_count'][0:len(exp[t]['N_count'])-1])
  nave[t] = np.mean(exp[t]['N_count'][0:len(exp[t]['N_count'])-1])

  covsum[t] = np.sum(exp[t]['cov_count'][0:len(exp[t]['cov_count'])-1])
  covave[t] = np.mean(exp[t]['cov_count'][0:len(exp[t]['cov_count'])-1])

  npsum[t] = np.sum(exp[t]['N_perc'][0:len(exp[t]['N_perc'])-1])
  npave[t] = np.mean(exp[t]['N_perc'][0:len(exp[t]['N_perc'])-1])

  cpsum[t] = np.sum(exp[t]['Cov_perc'][0:len(exp[t]['Cov_perc'])-1])
  cpave[t] = np.mean(exp[t]['Cov_perc'][0:len(exp[t]['Cov_perc'])-1])

for m in range(len(exp)):
  new_row = ['Total',l3sum[m],nsum[m],covsum[m],npsum[m],cpsum[m]]
  new_row2 = ['Average',l3ave[m],nave[m],covave[m],npave[m],cpave[m]]
  exp[m].loc["Total"] = new_row
  exp[m].loc["Average"] = new_row2
  ovr.append(cpave[m])
  tmp = list(exp[m]['Cov_perc'])
  del tmp[-3:]
  tmp2 = species.loc[m]['Species']
  gene[tmp2] = pd.Series(tmp)

gene = gene.set_index(['Gene'])
gene['Gene Mean'] = gene.mean(axis=1)
gene.loc['Species Mean'] = gene.mean()

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
