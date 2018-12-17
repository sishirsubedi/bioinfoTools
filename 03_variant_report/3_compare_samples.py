import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2

df_snp_1 = pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/COV_1_2_snp_indel_filter.tsv",sep='\t')
df_snp_mut_1 = df_snp_1[['chrom','position','ref','var']]

df_snp_2 = pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_1R_COV_2R/COV_1_2_R_snp_indel_filter.tsv",sep='\t')
df_snp_mut_2 = df_snp_2[['chrom','position','ref','var']]

df_join = pd.merge(df_snp_mut_1,df_snp_mut_2,on=['chrom','position','ref','var'],how='outer',indicator=True)

both = df_join[df_join['_merge']=='both'].shape[0]
sample_A = df_join[df_join['_merge']=='left_only'].shape[0]
sample_B = df_join[df_join['_merge']=='right_only'].shape[0]
venn2(subsets = (sample_A, sample_B, both), set_labels = ('COV_1_2', 'COV_1_2_R'))
plt.savefig('samples_venn')
plt.close()
