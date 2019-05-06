import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

# SAMPLE='COV-1_COV-2'
# SAMPLE='COV-1R_COV-2R'
SAMPLE='COLO-829_S5_COLO-829BL_S4'

###################### ADD refstd column to vep parsed file
refstd="/home/hhadmin/exome_pipeline/RefStd/reference_standard_creDesign.txt"
df_ref = pd.read_csv(refstd,sep='\t')

design= "/home/environments/ngs_test/exomeAnalysis/Paired/"+SAMPLE+"/"+SAMPLE+".variantcallers.combine.vep.parse.txt"
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


########################### generate count dictionary ######################
sample= "/home/environments/ngs_test/exomeAnalysis/Paired/"+SAMPLE+"/"+SAMPLE+".vc.combine.parsed_withREFSTD.vep.vcf.txt"
df_sample = pd.read_csv(sample,sep='\t')
df=df_sample
# df=df_sample[['VARSCAN','MUTECT', 'STRELKA', 'REFSTD']]
data={}
data['V']=str(df[ ( ( df['VARSCAN']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0])+';'\
         +str(df[df['VARSCAN']==1].shape[0] - df[ ( ( df['VARSCAN']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0] )
data['M']=str(df[ ( ( df['MUTECT']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0])+';'\
         +str(df[df['MUTECT']==1].shape[0] - df[ ( ( df['MUTECT']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0] )
data['S']=str(df[ ( ( df['STRELKA']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0])+';'\
          +str(df[df['STRELKA']==1].shape[0] - df[ ( ( df['STRELKA']==1 ) & ( df['REFSTD']==1 ) ) ].shape[0] )

## any one
df2=df[ df['REFSTD']==1]
anyone=df2[ ( df2['VARSCAN']==1 ) |  ( df2['MUTECT']==1 ) | ( df2['STRELKA']==1 )].shape[0]
data['X+R']=str(anyone) + ';' +str(df.shape[0] - anyone)


### all three

data['X+X+X+R']=str(df2[ ( df2['VARSCAN']==1 ) &  ( df2['MUTECT']==1 ) & ( df2['STRELKA']==1 )].shape[0]) + ';'\
            +str(df.shape[0] - df2[ ( df2['VARSCAN']==1 ) &  ( df2['MUTECT']==1 ) & ( df2['STRELKA']==1 )].shape[0] )


## any two
only_one=df[ ( \
                    ( df['VARSCAN']==1 ) &  ( df['MUTECT']==0 ) &  ( df['STRELKA']==0 )| \
                    ( df['MUTECT']==1 ) &  ( df['VARSCAN']==0 )&  ( df['STRELKA']==0 ) | \
                    ( df['STRELKA']==1 ) &  ( df['VARSCAN']==0 ) &  ( df['MUTECT']==0 )  \
                    )].index

df3=df.drop(only_one)

data['X+X+R']=str(df3[df3['REFSTD']==1].shape[0])+';'+str(df3.shape[0] - df3[df3['REFSTD']==1].shape[0] )


df3.to_csv("/home/environments/ngs_test/exomeAnalysis/Paired/"+SAMPLE+"/"+SAMPLE+".vc.combine.parsed_withREFSTD_Anytwo.vep.vcf.txt",sep='\t',index=False)



for val in data.keys():
    print(val,'\t',data[val])

########################################## generate vep result


sample= "/home/environments/ngs_test/exomeAnalysis/Paired/"+SAMPLE+"/"+SAMPLE+".vc.combine.parsed_withREFSTD_Anytwo.vep.vcf.txt"
df_sample = pd.read_csv(sample,sep='\t')

print(df_sample[df_sample['Consequence'].str.contains('missense_variant')].shape[0])
print(df_sample[df_sample['Consequence'].str.contains('inframe_insertion')].shape[0])
print(df_sample[df_sample['Consequence'].str.contains('inframe_deletion')].shape[0])
print(df_sample[df_sample['Consequence'].str.contains('stop_gained')].shape[0])
print(df_sample[df_sample['Consequence'].str.contains('frameshift_variant')].shape[0])
print(df_sample[df_sample['Consequence'].str.contains('coding_sequence_variant')].shape[0])
