import sys
import pandas as pd

file=sys.argv[1]
sample= file.split('.')[0]
print("processing..",sample)
df_snp = pd.read_excel(file,skiprows=5,sheetname=0)
# df_snp = df_snp[['Chromosome', 'Position (GRCh37)', 'Ref.', 'Alt.']]
df_snp = df_snp[['Chromosome', 'Position (GRCh37)', 'Ref.', 'Alt.']]
# df_snp.columns = ['CHROM','POS','REF','ALT']
df_snp.columns = ['CHROM','POS']
print(df_snp.shape)

df_indel = pd.read_excel(file,skiprows=5,sheetname=1)
# df_indel = df_indel[['Chromosome', 'Position (GRCh37)', 'Ref.', 'Alt.']]
df_indel = df_indel[['Chromosome', 'Position (GRCh37)', 'Ref.', 'Alt.']]
# df_indel.columns = ['CHROM','POS','REF','ALT']
df_indel.columns = ['CHROM','POS']
print(df_indel.shape)

df_combine = pd.concat([df_snp,df_indel])
print(df_combine.shape)
df_combine.drop_duplicates(keep='first',inplace=True)
print(df_combine.shape)
df_combine.to_csv(sample+".txt",sep='\t',index=False)


# combine all three
df1=pd.read_csv("2016-01-13_RKO_v3.txt",sep='\t')
print(df1.shape)
df2=pd.read_csv("2016-01-13_SW48_v3.txt",sep='\t')
print(df2.shape)
df3=pd.read_csv("2016-01-13_HCT116_v3.txt",sep='\t')
print(df3.shape)


df_combine = pd.concat([df1,df2,df3])
print(df_combine.shape)
df_combine.drop_duplicates(keep='first',inplace=True)
print(df_combine.shape)
df_combine.to_csv("R_S_H_Combine.txt",sep='\t',index=False)

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2

sample="Horizon"
sample1="hct116"
sample2="qm"

df_back = pd.read_csv("2016-01-13_HCT116_v3.txt",sep='\t')
df_back.head()
df_cell = pd.read_csv("QM_v3.txt",sep='\t')
df_cell.head()


df_join = pd.merge(df_back,df_cell,on=['CHROM','POS','REF','ALT'],how='outer',indicator=True)
both = df_join[df_join['_merge']=='both'].shape[0]
sample_A = df_join[df_join['_merge']=='left_only'].shape[0]
sample_B = df_join[df_join['_merge']=='right_only'].shape[0]
venn2(subsets = (sample_A, sample_B, both), set_labels = ( sample1,sample2))
plt.title(sample+" "+sample1+" "+sample2)
plt.savefig("venn_"+sample+"_"+sample1+"_"+sample2)
plt.close()

df_join.to_csv("horizon_comparison.txt",sep='\t',index=False)


#### confirm mutation

df_vmut = pd.read_csv("horizon_verified_mutations.txt",sep='\t')
# df_file= pd.read_csv("R_S_H_Combine.txt",sep='\t')
df_file= pd.read_csv("QM_v3.txt",sep='\t')


for i,r in df_vmut.iterrows():
    df_min = df_file[( (df_file['CHROM']==r['CHROM']) & (df_file['POS']==r['POS']) & (df_file['REF']==r['REF']) & (df_file['ALT']==r['ALT'])) ]
    if df_min.shape[0]>=1:
        print(df_min)
        print(r)
