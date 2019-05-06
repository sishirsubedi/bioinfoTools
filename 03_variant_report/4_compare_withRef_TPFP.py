import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib_venn import venn3

sample="COV1R"
sample1="varscan"
sample2="mutect"
sample3="strelka"

df_1 = pd.read_csv("/home/environments/ngs_test/exomeAnalysis/Paired/COLO-829_S5_COLO-829BL_S4/varscan/COLO-829_S5_COLO-829BL_S4.varscan.parafilter.crefilter",sep='\t')
df_mut_1 = df_1[['CHROM','POS','REF','ALT']]

df_2= pd.read_csv("/home/environments/ngs_test/exomeAnalysis/Paired/COLO-829_S5_COLO-829BL_S4/mutect/COLO-829_S5_COLO-829BL_S4.mutect.parafilter.crefilter",sep='\t')
df_mut_2 = df_2[['CHROM','POS','REF','ALT']]

df_3 = pd.read_csv("/home/environments/ngs_test/exomeAnalysis/Paired/COLO-829_S5_COLO-829BL_S4/strelka/COLO-829_S5_COLO-829BL_S4.strelka.parafilter.crefilter",sep='\t')
df_mut_3 = df_3[['CHROM','POS','REF','ALT']]


df_ref = pd.read_csv("/home/hhadmin/exome_pipeline/RefStd/reference_standard_creDesign.txt",sep='\t')
df_ref = df_ref[['CHROM','POS','REF','ALT']]

cols=[]
for indx,row in df_ref.iterrows():
    temp=[]
    if(df_mut_1[ ( (df_mut_1['CHROM']==row['CHROM']) & (df_mut_1['POS']==row['POS']) & (df_mut_1['REF']==row['REF']) &  (df_mut_1['ALT']==row['ALT']))].shape[0] == 1):
        temp.append(1)
    else:
        temp.append(0)

    if(df_mut_2[ ( (df_mut_2['CHROM']==row['CHROM']) & (df_mut_2['POS']==row['POS']) & (df_mut_2['REF']==row['REF']) &  (df_mut_2['ALT']==row['ALT']))].shape[0] == 1):
        temp.append(1)
    else:
        temp.append(0)

    if(df_mut_3[ ( (df_mut_3['CHROM']==row['CHROM']) & (df_mut_3['POS']==row['POS']) & (df_mut_3['REF']==row['REF']) &  (df_mut_3['ALT']==row['ALT']))].shape[0] == 1):
        temp.append(1)
    else:
        temp.append(0)


    cols.append(temp)

df_ref['varscan']= [x[0] for x in cols]
df_ref['mutect']= [x[1] for x in cols]
df_ref['strelka']= [x[2] for x in cols]

varscan = [df_ref['varscan'].sum(),df_mut_1.shape[0]-df_ref['varscan'].sum()]
mutect = [df_ref['mutect'].sum(),df_mut_2.shape[0]-df_ref['mutect'].sum()]
strelka = [df_ref['strelka'].sum(),df_mut_3.shape[0]-df_ref['strelka'].sum()]

barWidth = 1
r = [0,1,2]
# Create brown bars
bars1=[varscan[0],mutect[0],strelka[0]]
bars2=[varscan[1],mutect[1],strelka[1]]
plt.bar(r, bars1, color='#557f2d', edgecolor='white', width=barWidth,label='reference')
# Create green bars (middle), on top of the firs ones
plt.bar(r, bars2, bottom=bars1, color='#7f6d5f', edgecolor='white', width=barWidth,label='unique')
plt.xticks(r, ['varscan','mutect','strelka'], fontweight='bold')
plt.legend()
plt.xlabel("variant caller")
plt.ylabel("snps + indels")
plt.savefig("/home/environments/ngs_test/exomeAnalysis/Paired/COLO-829_S5_COLO-829BL_S4/"+sample+'_'+sample1+'_'+'_'+sample2+'_'+'_'+sample3+'_3bARPLOT')
plt.close()
