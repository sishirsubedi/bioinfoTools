import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2

sample="COV_1_COV_2"
sample1="ref_design"
sample2="varscan"
fig_name="venn_2_"+sample+"_"+sample1+"_"+sample2+"pval1percent"
# df_snp_1 = pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/reference_standard_creExonOnly.txt",sep='\t')
df_snp_1 = pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/reference_standard_creDesign.txt",sep='\t')
# df_snp_1 = pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/varscan/01_COV_1_COV_2_snp_indel_filter.vcf.txt",sep='\t')
# df_snp_1= pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/gatk/01_gatk_filter_combine.vcf.txt",sep='\t')
df_snp_mut_1 = df_snp_1[['CHROM','POS','REF','ALT']]

df_snp_2 =  pd.read_csv("/home/environments/ngs_test/exomeAnalysis/COV-1_COV-2/varscan/COV-1_COV-2_varscan_output.snp.indel.filter.1percent_th.vcf.txt",sep='\t')
# df_snp_2 = pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/gatk/01_gatk_filter_combine.vcf.txt",sep='\t')
# df_snp_2 = pd.read_csv("/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/strelka/results/variants/snvs.strelka.filter.txt",sep='\t')
df_snp_mut_2 = df_snp_2[['CHROM','POS','REF','ALT']]
print(df_snp_mut_2.shape)

df_join = pd.merge(left=df_snp_mut_1, right=df_snp_mut_2, on=['CHROM','POS','REF','ALT'],how='outer',indicator=True)
# df_join.to_csv("temp.txt",sep='\t')
both = df_join[df_join['_merge']=='both'].shape[0]
sample_A = df_join[df_join['_merge']=='left_only'].shape[0]
sample_B = df_join[df_join['_merge']=='right_only'].shape[0]
venn2(subsets = (sample_A, sample_B, both), set_labels = ( sample1,sample2))
plt.title(sample+" "+sample1+" "+sample2)
plt.savefig(fig_name)
plt.close()
