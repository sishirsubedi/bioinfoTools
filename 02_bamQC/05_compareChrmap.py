####### input arguments are:
## file1_chrmap file2_chrmap

### example
# /opt/python3/bin/python3 04_compareChrmap.py in_CRE_corriell/cre_chrmap

#outputs are two graphs named as:
# file1_file2_chr_map

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


df_file_1 = pd.read_csv(sys.argv[1],header=None,sep='\t')
print(df_file_1.shape)



df_file_1.columns = ['chromosome','chr_length','mapped','unmapped']
df_file_1['map_percent'] = df_file_1['mapped']/df_file_1['mapped'].sum() * 100
df_file_1['unmap_percent'] = df_file_1['unmapped']/df_file_1['unmapped'].sum() * 100
tagline_1=str(sys.argv[1]).split('/')
tag_1=tagline_1[len(tagline_1)-1].split('_')[0]
df_file_1['tag'] = tag_1


df_main = df_file_1

plt.rcParams["figure.figsize"] = (20,10)
plt.rcParams["font.size"] = (12)
plt.ylim(0, df_file_1['map_percent'].max()+2.5)
sns.barplot(x='chromosome',y='map_percent',hue='tag',data=df_main)
plt.savefig(tag_1+'_chr_map')
plt.close()


print('completed')
# sns.set(font_scale=2) # to zoom figure
# df_file_1.plot.bar(x='flag',y='count',ylim=(0,1000),color='red')
# df_file_2.plot.bar(x='flag',y='count',ylim=(0,1000),color='green')
