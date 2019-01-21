import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib.pyplot import figure

sample="COV_1_2"

with open(sample+"_snp_indel_filter.vep.vcf.txt") as myfile:
    head = [next(myfile) for x in range(2)]
line_2=str(head[1].strip().split(':')[1]).split('|')

df_vep = pd.read_csv(sample+"_snp_indel_filter.vep.vcf.txt",sep='\t',skiprows=2)
df_vep.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
df_vep[line_2] = df_vep['INFO'].str.split('|',expand=True)
df_vep.drop(['INFO'],axis=1,inplace=True)

df_vep= df_vep[df_vep['GIVEN_REF']==df_vep['USED_REF']]
df_vep.to_csv(sample+"_vep_output.txt",sep='\t',index=False)
