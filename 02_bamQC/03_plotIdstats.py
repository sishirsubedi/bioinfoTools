import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


df_idxstats = pd.read_csv(sys.argv[1],header=None,sep='\t')
df_idxstats.columns = ['chromosome','chr_length','mapped','unmapped']
filter_indx = df_idxstats[df_idxstats['chromosome'].str.contains("_")].index
df_idxstats.drop(filter_indx, inplace=True)
df_idxstats['map_percent'] = df_idxstats['mapped']/df_idxstats['mapped'].sum() * 100

tagline_1=str(sys.argv[1]).split('/')
tag_1=tagline_1[len(tagline_1)-1].split('.')[0]
df_idxstats['tag'] = tag_1

plt.rcParams["figure.figsize"] = (20,10)
plt.rcParams["font.size"] = (12)
plt.ylim(0, df_idxstats['map_percent'].max()+2.5)
sns.barplot(x='chromosome',y='map_percent',hue='tag',data=df_idxstats)
plt.ylabel('Percent of mapped sequences')
plt.savefig(tag_1+'.sorted.rmdups.filter.bam.idxstats.plot.png')
plt.close()
print('completed')


df_design = pd.read_csv(sys.argv[2],header=None,sep='\t')
df_design.columns = ['chromosome','start','end']
total = df_design.shape[0]
df_design['count'] = df_design.groupby('chromosome')['chromosome'].transform('count')
df_design.drop_duplicates('chromosome',keep='first',inplace=True)
df_design['count'] = df_design['count']/total * 100
plt.rcParams["figure.figsize"] = (20,10)
plt.rcParams["font.size"] = (12)
sns.barplot(x='chromosome',y='count',hue='chromosome',data=df_design)
plt.savefig('cre_design_coverage.bed.png')
plt.close()
