import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

SAMPLE='COLO-829_S5_COLO-829BL_S4'

###### ADD refstd column to vep parsed file
refstd="/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/reference_standard_creDesign.txt"
df_ref = pd.read_csv(refstd,sep='\t')

design= "/home/environments/ngs_test/exomeAnalysis/Paired/"+SAMPLE+"/"+SAMPLE+".vc.combine.parsed.vep.vcf.txt"
df_sample = pd.read_csv(design,sep='\t')
# df_sample.columns=['CHROM','START','END']

filter=[]
for indx,row in df_sample.iterrows():
    if ( df_ref[ (df_ref['CHROM'] == row['CHROM']) & \
        (df_ref['POS'] == row['POS']) & \
        (df_ref['REF'] == row['REF']) & \
        (df_ref['ALT'] == row['ALT']) ].shape[0]>=1):
        # print([row['CHROM'],row['POS'],row['REF'],row['ALT']])
        filter.append(1)
    else:
        filter.append(0)

df_sample['REFSTD']=filter
df_sample.to_csv("/home/environments/ngs_test/exomeAnalysis/Paired/"+SAMPLE+"/"+SAMPLE+".vc.combine.parsed_withREFSTD.vep.vcf.txt",sep='\t',index=False)


############################ generate count dictionary ######################
sample= "/home/environments/ngs_test/exomeAnalysis/Paired/"+SAMPLE+"/"+SAMPLE+".vc.combine.parsed_withREFSTD.vep.vcf.txt"
df_sample = pd.read_csv(sample,sep='\t')
df=df_sample[['VARSCAN','MUTECT', 'STRELKA', 'REFSTD']]
data={}
data['V']=df[df['VARSCAN']==1].shape[0]
data['M']=df[df['MUTECT']==1].shape[0]
data['S']=df[df['STRELKA']==1].shape[0]
data['R']=df[df['REFSTD']==1].shape[0]

data['V+R']=df[ ( ( df['VARSCAN']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0]
data['M+R']=df[ ( ( df['MUTECT']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0]
data['S+R']=df[ ( ( df['STRELKA']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0]

data['V+M']=df[ ( ( df['VARSCAN']==1 ) & ( df['MUTECT']==1 ) ) ].shape[0]
data['V+S']=df[ ( ( df['VARSCAN']==1 ) & ( df['STRELKA']==1 ) ) ].shape[0]
data['M+S']=df[ ( ( df['MUTECT']==1 ) & ( df['STRELKA']==1 ) ) ].shape[0]

data['V+M+R']=df[ ( ( df['VARSCAN']==1 ) &  ( df['MUTECT']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0]
data['V+S+R']=df[ ( ( df['VARSCAN']==1 ) &  ( df['STRELKA']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0]
data['M+S+R']=df[ ( ( df['MUTECT']==1 ) &  ( df['STRELKA']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0]
data['V+M+S']=df[ ( ( df['MUTECT']==1 ) &  ( df['STRELKA']==1 ) & ( df['VARSCAN']==1 ) ) ].shape[0]

data['V+M+S+R']=df[ ( ( df['VARSCAN']==1 ) &  ( df['MUTECT']==1 ) & ( df['STRELKA']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0]

df2=df[ df['REFSTD']==1]
data['X+R']=df2[ ( ( df2['VARSCAN']==1 ) |  ( df2['MUTECT']==1 ) | ( df2['STRELKA']==1 ))].shape[0]
data['X+X+R']=df2[ ( \
       ( df2['VARSCAN']==1 ) &  ( df2['MUTECT']==1 ) | \
       ( df2['VARSCAN']==1 ) &  ( df2['STRELKA']==1 ) | \
       ( df2['MUTECT']==1 ) &  ( df2['STRELKA']==1 )  \
       )].shape[0]

for val in data.keys():
    print(val,'\t',data[val])
