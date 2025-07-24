# RAY Script mosdepth
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

for i in range(len(names_foo)):
  h = pd.read_csv(names_foo[i], sep = '\t')
  exp.append(h)

lensum = np.zeros(len(exp))
basum = np.zeros(len(exp))
meansum = np.zeros(len(exp))
maxsum = np.zeros(len(exp))
minsum = np.zeros(len(exp))
lenave = np.zeros(len(exp))
basave = np.zeros(len(exp))
meanave = np.zeros(len(exp))
minave = np.zeros(len(exp))
maxave = np.zeros(len(exp))

for t in range(len(exp)):
  lensum[t] = np.sum(exp[t]['length'][0:len(exp[t]['length'])-1])
  lenave[t] = np.mean(exp[t]['length'][0:len(exp[t]['length'])-1])

  basum[t] = np.sum(exp[t]['bases'][0:len(exp[t]['bases'])-1])
  basave[t] = np.mean(exp[t]['bases'][0:len(exp[t]['bases'])-1])

  meansum[t] = np.sum(exp[t]['mean'][0:len(exp[t]['mean'])-1])
  meanave[t] = np.mean(exp[t]['mean'][0:len(exp[t]['mean'])-1])

  maxsum[t] = np.sum(exp[t]['max'][0:len(exp[t]['max'])-1])
  maxave[t] = np.mean(exp[t]['max'][0:len(exp[t]['max'])-1])

  minsum[t] = np.sum(exp[t]['min'][0:len(exp[t]['min'])-1])
  minave[t] = np.mean(exp[t]['min'][0:len(exp[t]['min'])-1])

for m in range(len(exp)):
  new_row = ['Total',lensum[m],basum[m],meansum[m],minsum[m],maxsum[m]]
  new_row2 = ['Average',lenave[m],basave[m],meanave[m],minave[m],maxave[m]]
  exp[m].loc["Total"] = new_row
  exp[m].loc["Average"] = new_row2

sheetname = []
for m in range(len(exp)):
  sheetname.append(df[0][m])

c = list(np.zeros(len(exp)))
for s in range(len(exp)):
  c[s] = df[0][s]

for i in range(len(exp)):
    exp[i].to_csv("%s.csv" % c[i], index=False)