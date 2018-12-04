import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib.pyplot import figure
figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')

df_file_1 = pd.read_csv(sys.argv[1],header=None,sep='\t')
df_file_1.columns = ['chromosome','design_start','design_end','read_depth']
tagline_1=str(sys.argv[1]).split('/')
tag_1=tagline_1[len(tagline_1)-1].split('.')[0]
df_file_1['tag'] = tag_1

plt.xlim(0,1000)
sns.distplot(df_file_1['read_depth'],bins=400)
plt.title(tag_1+'_Target coverage distribution')
plt.savefig(tag_1+'_uniformity_1')
plt.close()


figure(num=None, figsize=(20, 10), dpi=100, facecolor='w', edgecolor='k')
plt.ylim(0,1000)
plt.plot(df_file_1.index, df_file_1['read_depth'],'bo', markersize=0.5)
plt.xlabel("cre_targets")
plt.ylabel("read_depth")
plt.savefig(tag_1+'_uniformity_2')
plt.close()

def smooth(y, box_pts):
    box = np.ones(box_pts)*box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


figure(num=None, figsize=(20, 10), dpi=100, facecolor='w', edgecolor='k')
rolling = df_file_1['read_depth'].rolling(window=25)
rolling_mean = rolling.mean()
plt.ylim(0,1000)
plt.plot(df_file_1.index, rolling_mean,'g-', lw=1.0)
rolling = df_file_1['read_depth'].rolling(window=2500)
rolling_mean = rolling.mean()
plt.plot(df_file_1.index, rolling_mean,'r-', lw=1.0)
plt.xlabel("cre_targets")
plt.ylabel("smoothen_read_depth")
plt.savefig(tag_1+'_uniformity_3')
plt.close()

figure(num=None, figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')
df_file_10 = df_file_1[df_file_1['read_depth']<10]
plt.plot(range(df_file_10.shape[0]), df_file_10['read_depth'],'bo', markersize=2.0)
# plt.plot(range(df_file_10.shape[0]), smooth(df_file_10['read_depth'],1), 'g-', lw=1)
plt.xlabel("cre_targets")
plt.ylabel("read_depth")
plt.savefig(tag_1+'_uniformity_4')
plt.close()
