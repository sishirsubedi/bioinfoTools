import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib_venn import venn3,venn3_unweighted,venn3_circles

sample="COV"
sample1="strelka_COV1"
sample2="strelka_COV1R"
sample3="strelka_COLO"

df_1 = pd.read_csv("/home/environments/ngs_test/exomeAnalysis/Paired/COV-1_COV-2/strelka/COV-1_COV-2.strelka.filter.vcf.txt",sep='\t')
df_mut_1 = df_1[['CHROM','POS','REF','ALT']]

df_2= pd.read_csv("/home/environments/ngs_test/exomeAnalysis/Paired/COV-1R_COV-2R/strelka/COV-1R_COV-2R.strelka.filter.vcf.txt",sep='\t')
df_mut_2 = df_2[['CHROM','POS','REF','ALT']]

df_3 = pd.read_csv("/home/environments/ngs_test/exomeAnalysis/Paired/COLO-829_S5_COLO-829BL_S4/strelka/COLO-829_S5_COLO-829BL_S4.strelka.filter.vcf.txt",sep='\t')
df_mut_3 = df_3[['CHROM','POS','REF','ALT']]

df_join_1 = pd.merge(df_mut_1,df_mut_2,on=['CHROM','POS','REF','ALT'],how='outer',indicator=False)
df_join = pd.merge(df_join_1,df_mut_3,on=['CHROM','POS','REF','ALT'],how='outer',indicator=False)

print(df_join.shape)
cols=[]
for indx,row in df_join.iterrows():
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

df_join[sample1]= [x[0] for x in cols]
df_join[sample2]= [x[1] for x in cols]
df_join[sample3]= [x[2] for x in cols]

# df_join.to_csv("join_table.txt",sep='\t',index=False)

total = df_join.shape[0]
all_3 = df_join[(df_join[sample1]==1) & (df_join[sample2]==1) & (df_join[sample3]==1)].shape[0]
all_none = df_join[(df_join[sample1]==0) & (df_join[sample2]==0) & (df_join[sample3]==0)].shape[0]
all_2_cosmic_clinvar = df_join[(df_join[sample1]==1) & (df_join[sample2]==1) & (df_join[sample3]==0)].shape[0]
all_2_cosmic_gnomad = df_join[(df_join[sample1]==1) & (df_join[sample2]==0) & (df_join[sample3]==1)].shape[0]
all_2_clinvar_gnomad =df_join[(df_join[sample1]==0) & (df_join[sample2]==1) & (df_join[sample3]==1)].shape[0]
unq_cosmic = df_join[(df_join[sample1]==1) & (df_join[sample2]==0) & (df_join[sample3]==0)].shape[0]
unq_clinvar = df_join[(df_join[sample1]==0) & (df_join[sample2]==1) & (df_join[sample3]==0)].shape[0]
unq_gnomad = df_join[(df_join[sample1]==0) & (df_join[sample2]==0) & (df_join[sample3]==1)].shape[0]
unknown = total - all_3- all_2_cosmic_clinvar - all_2_cosmic_gnomad - all_2_clinvar_gnomad - unq_cosmic - unq_clinvar - unq_gnomad

plt.figure(num=None, figsize=(10, 8), dpi=80, facecolor='w', edgecolor='k')
out = venn3(subsets = (unq_cosmic,unq_clinvar, all_2_cosmic_clinvar, unq_gnomad,all_2_cosmic_gnomad,all_2_clinvar_gnomad,all_3), \
    set_labels = (sample1, sample2, sample3))
# for text in out.set_labels:
#     text.set_fontsize(20)
# for text in out.subset_labels:
#     text.set_fontsize(12)

    # id2idx = {'10': 0, '01': 1, '11': 2,
              # '100': 0, '010': 1, '110': 2, '001': 3, '101': 4, '011': 5, '111': 6, 'A': 0, 'B': 1, 'C': 2}
# label = out.get_label_by_id('11')          # Those are subset labels (i.e. numbers)
# label.set_x(label.get_position()[0] + 0.1)
plt.savefig(sample+'_'+sample1+'_'+'_'+sample2+'_'+'_'+sample3+'_3venn')
plt.close()
