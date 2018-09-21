####### input arguments are:
## file1_flagstat_mapq file2_flagstat_mapq graph
#here graph can be 'flag' or 'mapq'

### example
# /opt/python3/bin/python3 04_compareMapq.py in_CRE_corriell/cre_samflag_mapq in_V7_corriell/XTHS_100ng_hg19/v7_samflag_mapq flag

#outputs are two graphs named as:
# graph_1percent and graph_100percent

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


df_file_1 = pd.read_csv(sys.argv[1],header=None,sep='\t')
compare = sys.argv[2]

print(df_file_1.shape)
print(compare)

df_file_1.columns = ['flag','mapq']
df_file_1['count'] = df_file_1.groupby(compare)[compare].transform('count')
df_file_1['count'] = df_file_1['count']/df_file_1.shape[0] * 100
df_file_1.drop_duplicates(compare,keep='first',inplace=True)
df_file_1.reset_index(inplace=True)
df_file_1.drop(['index'], 1,inplace=True)
df_file_1.sort_values(compare,inplace=True)
tagline_1=str(sys.argv[1]).split('/')
tag_1=tagline_1[len(tagline_1)-1].split('_')[0]
df_file_1['tag'] = tag_1


df_main = df_file_1

plt.rcParams["figure.figsize"] = (20,10)
plt.rcParams["font.size"] = (12)
plt.ylim(0, 100)
sns.barplot(x=compare,y='count',hue='tag',data=df_main)
plt.savefig(tag_1 +'_'+compare+'_100percent2')
plt.close()

plt.rcParams["figure.figsize"] = (20,10)
plt.rcParams["font.size"] = (12)
plt.ylim(0, 1)
sns.barplot(x=compare,y='count',hue='tag',data=df_main)
plt.savefig(tag_1 +'_'+compare+'_1percent2')
plt.close()


print('completed')
# sns.set(font_scale=2) # to zoom figure
# df_file_1.plot.bar(x='flag',y='count',ylim=(0,1000),color='red')
# df_file_2.plot.bar(x='flag',y='count',ylim=(0,1000),color='green')
