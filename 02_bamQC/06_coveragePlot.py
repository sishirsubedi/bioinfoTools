####### input arguments are:
## cre_bedcoverage

### example
# /opt/python3/bin/python3 04_compareChrmap.py in_CRE_corriell/cre_bedcoverage

#outputs are two graphs named as:
# file1_chromosome_coverage_plot


import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


df_file_1 = pd.read_csv(sys.argv[1],header=None,sep='\t')
chromosome = sys.argv[2]

print(df_file_1.shape)
print(chromosome)

df_file_1.columns = ['chromosome','design_start','design_end','design_base_num','count']
df_file_1 = df_file_1[df_file_1['chromosome']==chromosome]
tagline_1=str(sys.argv[1]).split('/')
tag_1=tagline_1[len(tagline_1)-1].split('_')[0]
df_file_1['tag'] = tag_1

df_main= df_file_1
cord = chromosome+'_cordinates'
df_main[cord] = 0
df_main[cord][df_main['tag']==tag_1] = list(range(len(df_main[cord][df_main['tag']==tag_1])))

df_main = df_main.iloc[1:1000,:]

plt.rcParams["figure.figsize"] = (20,10)
plt.rcParams["font.size"] = (12)
plt.ylim(0, df_main['count'].max()+2.5)
sns.lineplot(x=cord,y='count',hue='tag',data=df_main)
plt.savefig(tag_1+'_'+chromosome+'_coverage_plot')
plt.close()
print('completed')
