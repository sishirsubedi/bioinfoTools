####### input arguments are:
## cre_bedcoverage_mean

### example
# /opt/python3/bin/python3 04_compareChrmap.py in_CRE_corriell/cre_bedcoverage_mean

#outputs are two graphs named as:
# file1_x_coverage_QC

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np


df_file_1 = pd.read_csv(sys.argv[1],header=None,sep='\t')
coverage = int(sys.argv[2])

print(df_file_1.shape)
print('coverage is ' + str(coverage) + 'x')


df_file_1.columns = ['chromosome','design_start','design_end','count']
tagline_1=str(sys.argv[1]).split('/')
tag_1=tagline_1[len(tagline_1)-1].split('_')[0]
df_file_1['tag'] = tag_1

result=[]
covx = np.arange(1,coverage, 1)
for x in covx:
    result.append([x,(sum(df_file_1['count']>x)/df_file_1.shape[0]) * 100,tag_1])

df_main = pd.DataFrame(result)
df_main.columns = ['read_depth','average_count','tag']
plt.rcParams["figure.figsize"] = (20,10)
plt.rcParams["font.size"] = (12)
plt.ylim(0,  max(df_main['average_count'].max(),df_main['average_count'].max())+2.5)
plt.xlim(1,coverage)
sns.lineplot(x='read_depth',y='average_count',hue='tag',data=df_main)
plt.xlabel("Cumulative Read Depth")
plt.ylabel("% of Covered Target bases ")
plt.savefig(tag_1+'_'+str(coverage)+'x_coverage_QC')
plt.close()
print('completed')
