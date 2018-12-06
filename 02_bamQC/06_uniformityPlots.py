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

df_file_b10 = df_file_1[df_file_1['read_depth']<10]
df_file_b10.to_csv("COV-1.CREcoverage.mean.b10x.bed",sep='\t',header=False,index=False)

df_a10 = df_file_1[(df_file_1['read_depth']>100) & (df_file_1['read_depth']<300) ]
df_a10_sample = df_a10.sample(n=12114, replace=False)
df_a10_sample.to_csv("COV-1.CREcoverage.mean.above10x.bed",sep='\t',header=False,index=False)

/opt/bedtools2/bin/bedtools getfasta -fi /home/doc/ref/ref_genome/ucsc.hg19.fasta -bed COV-1.CREcoverage.mean.b10x.bed -tab -fo COV-1.CREcoverage.mean.b10x.fa.out
/opt/bedtools2/bin/bedtools getfasta -fi /home/doc/ref/ref_genome/ucsc.hg19.fasta -bed COV-1.CREcoverage.mean.above10x.bed -tab -fo COV-1.CREcoverage.mean.above10x.fa.out


df_b10 = pd.read_csv("COV-1.CREcoverage.mean.b10x.fa.out",sep='\t',header=None)
df_b10.columns=['region','seq']
at_count=[]
for indx,row in df_b10.iterrows():
    A = row['seq'].lower().count('a')
    T = row['seq'].lower().count('t')
    G = row['seq'].lower().count('g')
    C = row['seq'].lower().count('c')
    at_count.append(((A+T)/(A+T+G+C)) * 100)
df_b10['at_count']=at_count


df_above10 = pd.read_csv("COV-1.CREcoverage.mean.above10x.fa.out",sep='\t',header=None)
df_above10.columns=['region','seq']
at_count=[]
for indx,row in df_above10.iterrows():
    A = row['seq'].lower().count('a')
    T = row['seq'].lower().count('t')
    G = row['seq'].lower().count('g')
    C = row['seq'].lower().count('c')
    at_count.append(((A+T)/(A+T+G+C)) * 100)
df_above10['at_count']=at_count

th_mean=df_above10['at_count'].mean()
th_var=df_above10['at_count'].var()
th_samples = np.random.normal(th_mean, 6.0, 12114)

sns.distplot( df_b10['at_count'] , color="red", label="AT% < 10x",hist=False)
sns.distplot( df_above10['at_count'] , color="green", label="AT% (100~300)x",hist=False)
sns.distplot( th_samples , color="blue", label="AT% Theoretical",hist=False)
plt.title("AT distribution over sequences")
plt.xlabel("Mean AT %")
plt.ylabel("Proportion of reads")
plt.legend()
plt.savefig('ATrich_uniformity_test')
plt.close()
