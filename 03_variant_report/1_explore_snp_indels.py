import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

sample='COV_7_COV_8'
df_snp = pd.read_csv("output.snp.txt",sep='\t')
df_snp['chrom'] =df_snp['chrom'].apply(str)
df_snp['group']='snp'

df_indel = pd.read_csv("output.indel.txt",sep='\t')
df_indel['chrom'] =df_indel['chrom'].apply(str)
df_indel['group']='indel'
for indx, row in df_indel.iterrows():
    if row['var'][0] == '-':
        df_indel.set_value(indx, "ref", row['ref']+row['var'][1:])
        df_indel.set_value(indx, "var", row['ref'])
    elif row['var'][0] == '+':
        df_indel.set_value(indx, "var", row['ref'] + row['var'][1:] )

df_both = pd.concat([df_snp,df_indel])


df_both.sort_values(by=['somatic_status'],inplace=True)
sns.countplot(x=df_both['somatic_status'],hue=df_both['group'],data=pd.melt(df_both))
plt.savefig('snp_indel_count_raw')
plt.close()
print('plot completed')

###filtering:
filter_indx = df_both[df_both['chrom'].str.contains("GL|gl")].index
df_both.drop(filter_indx, inplace=True)
print('Raw:'+ str(df_both.shape[0]))
print('Raw:snp:'+ str(df_both[df_both['group']=='snp'].shape[0]))
print('Raw:indel:'+ str(df_both[df_both['group']=='indel'].shape[0]))


# somatic only
df_both_pass_s0 = df_both[(df_both['somatic_status'].isin(['Somatic']) )]
print('Somatic:'+ str(df_both_pass_s0.shape[0]))
print('Somatic:snp:'+ str(df_both_pass_s0[df_both_pass_s0['group']=='snp'].shape[0]))
print('Somatc:indel:'+ str(df_both_pass_s0[df_both_pass_s0['group']=='indel'].shape[0]))

# pval
df_both_pass_s1 = df_both_pass_s0[ (df_both_pass_s0['somatic_p_value']<0.05)]

# normal var freq is zero
df_both_pass_s2 = df_both_pass_s1[df_both_pass_s1['normal_var_freq'].str.replace('%','' ).astype(float)<=0]

# min depth 20x for both normal and tumor_reads
df_both_pass_s2['normal_depth'] = df_both_pass_s2['normal_reads1'] + df_both_pass_s2['normal_reads2']
df_both_pass_s2['tumor_depth'] = df_both_pass_s2['tumor_reads1'] + df_both_pass_s2['tumor_reads2']
df_both_pass_s3 = df_both_pass_s2[ (df_both_pass_s2['normal_depth']>=20) & (df_both_pass_s2['tumor_depth']>=20) ]

print('Filter:'+ str(df_both_pass_s3.shape[0]))
print('Filter:snp:'+ str(df_both_pass_s3[df_both_pass_s3['group']=='snp'].shape[0]))
print('Filter:indel:'+ str(df_both_pass_s3[df_both_pass_s3['group']=='indel'].shape[0]))

vcf=[]
for indx,row in df_both_pass_s3.iterrows():
    info='ND:'+str(row['normal_depth'])+','+'TD:'+str(row['tumor_depth'])
    vcf.append([row['chrom'],row['position'],'.',row['ref'],row['var'],'.','.',info])

df_vcf=pd.DataFrame(vcf)
df_vcf.columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
df_vcf.to_csv("01_"+sample+"_snp_indel_filter.vcf.txt",sep='\t',index=False)
