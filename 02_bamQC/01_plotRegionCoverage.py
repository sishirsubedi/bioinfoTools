import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


df_bed = pd.read_csv(sys.argv[1],header=None,sep='\t')
df_bed.columns = ['chromosome','design_start','design_end','design_base_num','count']

df_genelist = pd.read_csv(sys.argv[2])

for indx,row in df_genelist.iterrows():
    df_bed_gene = df_bed[df_bed['chromosome']==chromosome]
    df_bed['tag'] = row['gene']

df_main= df_bed
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
