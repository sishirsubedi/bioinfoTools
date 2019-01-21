import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2

sample="COV_7_COV_8"
sample1="mutect"
sample2="varscan"
df_snp_1 = pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_7_COV_8/gatk/01_gatk_filter_combine.vcf.txt",sep='\t')
# df_snp_1= pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/gatk/01_gatk_filter_combine.vcf.txt",sep='\t')
df_snp_mut_1 = df_snp_1[['CHROM','POS','REF','ALT']]

df_snp_2 = pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_7_COV_8/01_COV_7_COV_8_snp_indel_filter.vcf.txt",sep='\t')
# df_snp_2 = pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/varscan/01_COV_1_COV_2_snp_indel_filter.vcf.txt",sep='\t')
df_snp_mut_2 = df_snp_2[['CHROM','POS','REF','ALT']]
print(df_snp_mut_2.shape)

df_join = pd.merge(df_snp_mut_1,df_snp_mut_2,on=['CHROM','POS','REF','ALT'],how='outer',indicator=True)
both = df_join[df_join['_merge']=='both'].shape[0]
sample_A = df_join[df_join['_merge']=='left_only'].shape[0]
sample_B = df_join[df_join['_merge']=='right_only'].shape[0]
venn2(subsets = (sample_A, sample_B, both), set_labels = ( sample1,sample2))
plt.title(sample+" "+sample1+" "+sample2)
plt.savefig("venn_"+sample+"_"+sample1+"_"+sample2)
plt.close()
