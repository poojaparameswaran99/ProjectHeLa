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


AnnotatedKinaseSheet = f'../10415/10415_unbiased_normalized_041124_AnnotatedHumanKinome.csv'
PhosphoSheet = f'../10415/10415_phos_normalized_041124.xlsx'

data = pd.read_csv(AnnotatedKinaseSheet, index_col=0)
phospho = pd.read_excel(f'{PhosphoSheet}', sheet_name='Normalized Data')

data['Human_Kinase'] = data['Human_Kinase'].replace('False', False)

kinases = data[data['Human_Kinase'] != False].reset_index(drop=True)

out = list(set(phospho['Accession']) & set(kinases['Accession']))
print(len(out), kinases.shape)

# make non-phospho kinases unique
if '_k' not in kinases.loc[0, 'Accession'].lower():
    kinases['Accession'] = kinases['Accession'].apply(lambda x: f'{x.split(";")[0]}_k')

all = pd.concat([kinases, phospho], axis=0)
print(f'Phospho shape: {phospho.shape}\nKinase shape {kinases.shape}\nConcatenatedFinal: {all.shape}')
all.to_csv(f'{os.path.dirname(AnnotatedKinaseSheet)}/KinSub{os.path.basename(AnnotatedKinaseSheet)}')
