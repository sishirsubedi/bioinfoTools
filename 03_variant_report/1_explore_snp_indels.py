import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

df_snp = pd.read_csv("output.snp.txt",sep='\t')
df_snp['group']='snp'

df_indel = pd.read_csv("output.indel.txt",sep='\t')
df_indel['group']='indel'
for indx, row in df_indel.iterrows():
    if row['var'][0] == '-':
        df_indel.set_value(indx, "ref", row['ref']+row['var'][1:])
        df_indel.set_value(indx, "var", row['ref'])
    elif row['var'][0] == '+':
        df_indel.set_value(indx, "var", row['ref'] + row['var'][1:] )



df_both = pd.concat([df_snp,df_indel])
sns.countplot(x=df_both['somatic_status'],hue=df_both['group'],data=pd.melt(df_both))
plt.savefig('snp_indel_count_raw')
plt.close()
print('plot completed')

df_both_pass = df_both[ (df_both['somatic_p_value']<0.01) & df_both['somatic_status'].isin(['Somatic','LOH'])]
sns.countplot(x=df_both_pass['somatic_status'],hue=df_both_pass['group'],data=pd.melt(df_both_pass))
plt.savefig('snp_indel_count_pval')
plt.close()
print('plot completed')
