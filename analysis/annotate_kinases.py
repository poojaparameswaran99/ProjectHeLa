import pandas as pd 
import numpy as np 
import os
import re
import matplotlib.pyplot as plt
from Bio import SeqIO
from itertools import groupby
import seaborn as sns
currdir = os.getcwd()
parent = os.path.dirname(currdir)
gparent = os.path.dirname(parent)

def ensure_dirs_exists(path):
    if "." in path:
	    path = os.path.dirname(path)
    if not os.path.exists(path):
        os.makedirs(path)
    return

## read the module parittion file in data/
mfile = f"../10415/10415_unbiased_normalized_041124.xlsx"
data = pd.read_excel(mfile, sheet_name='Normalized Data')
print(data.head())

## map human kinases
kinomefile = f"/home/poojaparameswaran/Documents/generalData/data/pkinfam.csv"
allkinome = pd.read_csv(kinomefile, index_col=0)
## compress into triples. 
humans = {}
for i, r in allkinome.iterrows():
    humans[r['HumanID']] = (r['HumanName'], r['ProteinName'])
print(humans)
data.columns
todrop = [x for x in data.columns if '_CV' in x]
if any('_CV' in col for col in data.columns):
    data.drop(columns=todrop, inplace=True)
s1 = data.copy(deep=True)

if 'Human_Kinase' not in s1.columns:
    s1.insert(s1.columns.get_loc('Accession')+1, 'Human_Kinase', '')
for i, r in s1.iterrows():
    poss = r['Accession'].split(';')
    match = False
    for x in poss:
        if x in humans.keys():
            s1.at[i, 'Human_Kinase'] = humans[x][0]
            match=True
            break
    if not match:
        s1.at[i, 'Human_Kinase'] = False
s2 = s1[s1['Human_Kinase'] != False]
s1.to_csv(f'{mfile.split(".xlsx")[0]}_AnnotatedHumanKinome.csv')
