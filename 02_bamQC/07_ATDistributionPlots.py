import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pandas as pd
import seaborn as sns
import numpy as np
import subprocess

def normalize(df,gap):
    ll=0.0
    ul=0.0+gap
    temp=[]
    while ul<=1.0 :
        val=df[(df['AT_Percent']>=ll) & (df['AT_Percent']<ul)].shape[0]/df.shape[0]
        ll += gap
        ul += gap
        temp.append(val)
    return temp


SAMPLE_1_BED_ATcount = sys.argv[1]
SAMPLE_2_BED_ATcount = sys.argv[2]
OUT_DIR = sys.argv[3]
OUTPUT_DIR_SAMPLE=os.path.abspath(os.path.dirname(SAMPLE_2_BED_ATcount))
SAMPLE_1=os.path.basename(SAMPLE_1_BED_ATcount).split('.')[0]
SAMPLE_2=os.path.basename(SAMPLE_2_BED_ATcount).split('.')[0]
print(SAMPLE_1)
print(SAMPLE_2)
print(OUTPUT_DIR_SAMPLE)


df_sample_1 = pd.read_csv(SAMPLE_1_BED_ATcount,header=None,sep='\t')
df_sample_1.columns =['AT_Percent']
print(df_sample_1.head())


df_sample_2 = pd.read_csv(SAMPLE_2_BED_ATcount,header=None,sep='\t')
df_sample_2.columns =['AT_Percent']
print(df_sample_2.head())


gap=0.025
df_norm = pd.DataFrame()
df_norm[SAMPLE_1]=normalize(df_sample_1,gap)
df_norm[SAMPLE_2]=normalize(df_sample_2,gap)
df_norm.index = np.arange(39)


df_norm.plot.bar(figsize=(15,10))
plt.xlabel("AT Content (%)")
plt.ylabel("Sequence Proportion")
plt.savefig(OUT_DIR+SAMPLE_1+"_"+SAMPLE_2+"_AT_Distribution_Norm")
plt.close()
