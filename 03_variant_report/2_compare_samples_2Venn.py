import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2


SAMPLE_1_FILE = sys.argv[1]
SAMPLE_2_FILE = sys.argv[2]
METHODS = sys.argv[3]
OUTPUT_DIR_SAMPLE=os.path.abspath(os.path.dirname(SAMPLE_1_FILE))
SAMPLE_1=os.path.basename(SAMPLE_1_FILE).split('.')[0]
METHOD_1 = METHODS.split('-')[0]
SAMPLE_2=os.path.basename(SAMPLE_2_FILE).split('.')[0]
METHOD_2 = METHODS.split('-')[1]
print(SAMPLE_1, METHOD_1)
print(SAMPLE_2, METHOD_2)
print(OUTPUT_DIR_SAMPLE)

FIG_NAME = OUTPUT_DIR_SAMPLE+"/"+SAMPLE_1+"_"+SAMPLE_2+"_"+METHOD_1+"_"+METHOD_2+"_VENN_2"

df_snp_1 =  pd.read_csv(SAMPLE_1_FILE,sep='\t')
df_snp_mut_1 = df_snp_1[['CHROM','POS','REF','ALT']]

df_snp_2 = pd.read_csv(SAMPLE_2_FILE,sep='\t')
df_snp_mut_2 = df_snp_2[['CHROM','POS','REF','ALT']]


df_join = pd.merge(left=df_snp_mut_1, right=df_snp_mut_2, on=['CHROM','POS','REF','ALT'],how='outer',indicator=True)

df_join.to_csv(OUTPUT_DIR_SAMPLE+'/'+SAMPLE_1+".2VENN.OUT.txt",sep='\t',index=False)


BOTH = df_join[df_join['_merge']=='both'].shape[0]
METHOD_1_NUM = df_join[df_join['_merge']=='left_only'].shape[0]
METHOD_2_NUM = df_join[df_join['_merge']=='right_only'].shape[0]

venn2(subsets = ( METHOD_2_NUM,METHOD_1_NUM, BOTH), set_labels = ( SAMPLE_2+" "+METHOD_2,SAMPLE_1+" "+METHOD_1))
plt.title(SAMPLE_1+" -- "+SAMPLE_2)
plt.savefig(FIG_NAME)
plt.close()
