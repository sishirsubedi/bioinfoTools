import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


df_file_1 = pd.read_csv(sys.argv[1],header=None)
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
tagline_1=str(sys.argv[1]).split('.')
tag_1=tagline_1[0]
df_file_1['tag'] = tag_1

plt.rcParams["figure.figsize"] = (20,10)
plt.rcParams["font.size"] = (12)
plt.ylim(0, 100)
plt.ylabel('Percentage of Framgents')
sns.barplot(x=compare,y='count',hue='tag',data=df_file_1)
plt.savefig(tag_1 +'.sorted.rmdups.filter.bam'+compare+'.100percent.png')
plt.close()

plt.rcParams["figure.figsize"] = (20,10)
plt.rcParams["font.size"] = (12)
plt.ylim(0, 1)
plt.ylabel('Percentage of Framgents')
sns.barplot(x=compare,y='count',hue='tag',data=df_file_1)
plt.savefig(tag_1 +'.sorted.rmdups.filter.bam'+compare+'.1percent.png')
plt.close()

print('completed')
