
############ Coverage method ############################

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

def checkPosition(x1,y1,x2,y2,olap):
    x1 += olap
    y1 -= olap
    if (x2 >= x1 and x2 <= y1) or (y2 >= x1 and y2 <= y1) :
        return True
    else:
        return False

def filterExons(df_e,df_d,olap,condition):
    filter=[]
    for indx,row in df_d.iterrows():
        df_e_mini=pd.DataFrame()
        if condition=='pure':
            df_e_mini = df_e[(df_e['chr']==row['chr']) &
            ( ( row['start'] >= df_e['start']) & ( row['end'] <= df_e['end'] ) ) ]
        elif condition=='olap':
            df_e_mini = df_e[(df_e['chr']==row['chr']) &
            ( ( row['start'] >= (df_e['start'])) & ( row['start'] <= (df_e['end']-olap) ) ) |
            ( ( row['end'] >= (df_e['start']+olap)) & ( row['end'] <= (df_e['end']) ) ) ]

        if df_e_mini.shape[0]>=1:
            print([row['chr'],row['start'],row['end'],df_e_mini.shape[0]])
            filter.append([row['chr'],row['start'],row['end'],df_e_mini.shape[0]])
    return filter

exons= "/home/hhadmin/exome_pipeline/01_bamQC/ucsc_hg19_exons_only_coordinates.txt"
design= "/home/hhadmin/exome_pipeline/01_bamQC/cre_design.bed"
output="/home/hhadmin/exome_pipeline/01_bamQC/"

df_exons= pd.read_csv(exons,sep='\t',header=None)
df_exons=df_exons.iloc[:,0:3]
df_exons.columns=["chr","start",'end']
filter_indx = df_exons[df_exons['chr'].str.contains("_")].index
df_exons.drop(filter_indx, inplace=True)
df_exons.set_index(['chr','start','end'])
# df_exons['exon_len']= df_exons['end']-df_exons['start']
# sns.distplot(df_exons[df_exons['exon_len']<1000]['exon_len'])
# sns.distplot(df_exons['exon_len'])
# plt.savefig('exons_len2')
# plt.close()

df_design= pd.read_csv(design,sep='\t',header=None)
df_design.columns=["chr","start",'end']
df_design.set_index(['chr','start','end'])
# df_design['design_len']= df_design['end']-df_design['start']
# sns.distplot(df_design['design_len'])
# # sns.distplot(df_design[df_design['design_len']<1000]['design_len'])
# plt.savefig('design_len_2')
# plt.close()

result= filterExons(df_exons,df_design,150,'pure')
pd.DataFrame(result).to_csv(output+"cre_design_exonOnly_withexonCount_pure.bed.txt",sep='\t',index=False)



############ Overlap method ############################

## use bedtools to get overlab
# /opt/bedtools2/bin/bedtools intersect -a cre_design.bed  -b ucsc_hg19_exons_only_coordinates_3columns.txt -wo > crev2_design_ucsc_exon_intersect.txt

#use pandas to get min,max,average overlapS
import pandas as pd
df = pd.read_csv("crev2_design_ucsc_exon_intersect.txt",sep='\t',header=None)
df = df.iloc[:,[0,1,2,6]]
df.columns=['chr','start','end','coverage']
df.groupby(['chr','start','end'])
df2 = df.groupby(['chr','start','end'])
df2 = df.groupby(['chr','start','end']).max()
df2_max = df.groupby(['chr','start','end']).max()
df2_min = df.groupby(['chr','start','end']).min()
df2_average = df.groupby(['chr','start','end']).mean()
df2_max['coverage'].sum()
df2_min['coverage'].sum()
df2_average['coverage'].sum()
