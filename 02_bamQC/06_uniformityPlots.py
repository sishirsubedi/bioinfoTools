import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pandas as pd
import seaborn as sns
import numpy as np
import subprocess

def plotHistogram(df,tag):
    figure(num=None, figsize=(20, 10), dpi=100, facecolor='w', edgecolor='k')
    plt.xlim(0,1000)
    sns.distplot(df['read_depth'],bins=400)
    plt.title(tag+'_Target coverage distribution')
    plt.savefig(tag+'_uniformity_depth_histogram')
    plt.close()

def smooth(y, box_pts):
    box = np.ones(box_pts)*box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def plotRawDepth(df,tag,smooth):
    figure(num=None, figsize=(20, 10), dpi=100, facecolor='w', edgecolor='k')
    plt.ylim(0,df['read_depth'].max())
    plt.plot(df.index, df['read_depth'],'bo', markersize=0.5)
    plt.xlabel("cre_targets")
    plt.ylabel("read_depth")
    plt.savefig(tag+'_uniformity_rawdepth')
    plt.close()

    if smooth:
        figure(num=None, figsize=(20, 10), dpi=100, facecolor='w', edgecolor='k')
        rolling = df['read_depth'].rolling(window=25)
        rolling_mean = rolling.mean()
        plt.ylim(0,df['read_depth'].max())
        plt.plot(df.index, rolling_mean,'g-', lw=1.0)
        rolling = df['read_depth'].rolling(window=2500)
        rolling_mean = rolling.mean()
        plt.plot(df.index, rolling_mean,'r-', lw=1.0)
        plt.xlabel("cre_targets")
        plt.ylabel("smoothen_read_depth")
        plt.savefig(tag_1+'_uniformity_smooth_rawdepth')
        plt.close()


def executeCommand(input,output):
    command = "/opt/bedtools2/bin/bedtools getfasta -fi /home/doc/ref/ref_genome/ucsc.hg19.fasta -bed {} -tab -fo {} ".format(input,output)
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return error

def getATCount(df):
    df.columns=['region','seq']
    at_count=[]
    for indx,row in df.iterrows():
        A = row['seq'].lower().count('a')
        T = row['seq'].lower().count('t')
        G = row['seq'].lower().count('g')
        C = row['seq'].lower().count('c')
        at_count.append(((A+T)/(A+T+G+C)) * 100)
    df['at_count']=at_count
    return df


df_coverage_all = pd.read_csv(sys.argv[1],header=None,sep='\t')
df_coverage_all.columns = ['chromosome','design_start','design_end','read_depth']
tagline_1=str(sys.argv[1]).split('/')
tag_1=tagline_1[len(tagline_1)-1].split('.')[0]
df_coverage_all['tag'] = tag_1

plotHistogram(df_coverage_all,tag_1)
plotRawDepth(df_coverage_all,tag_1,1)

df_coverage_below_10x = df_coverage_all[df_coverage_all['read_depth']<10]
plotRawDepth(df_coverage_below_10x,tag_1+'_below_10x',0)
df_coverage_below_10x.to_csv(tag_1+".CREcoverage.mean.b10x.bed",sep='\t',header=False,index=False)

# df_coverage_above_10x = df_coverage_all[(df_coverage_all['read_depth']>100) & (df_coverage_all['read_depth']<270) ]
df_coverage_above_10x = df_coverage_all[(df_coverage_all['read_depth']>10)]
# df_a10_sample = df_coverage_above_10x.sample(n=df_coverage_below_10x.shape[0], replace=False)
df_a10_sample = df_coverage_above_10x
df_a10_sample.to_csv(tag_1+".CREcoverage.mean.above10x.bed",sep='\t',header=False,index=False)


df_coverage_all.to_csv(tag_1+".CREcoverage.mean.allx.bed",sep='\t',header=False,index=False)

###create fasta file for specific regions
executeCommand(tag_1+".CREcoverage.mean.b10x.bed", tag_1+".CREcoverage.mean.b10x.fa.out")
executeCommand(tag_1+".CREcoverage.mean.above10x.bed", tag_1+".CREcoverage.mean.above10x.fa.out")
executeCommand(tag_1+".CREcoverage.mean.allx.bed", tag_1+".CREcoverage.mean.allx.fa.out")

df_b10 = pd.read_csv(tag_1+".CREcoverage.mean.b10x.fa.out",sep='\t',header=None)
getATCount(df_b10)

df_above10 = pd.read_csv(tag_1+".CREcoverage.mean.above10x.fa.out",sep='\t',header=None)
getATCount(df_above10)

# th_mean=df_above10['at_count'].mean()
# th_var=df_above10['at_count'].var()
# th_samples = np.random.normal(th_mean, 6.0, 12114)
#
# plt.xlim(0,100)
# sns.distplot( df_b10['at_count'] , color="red", label="AT% < 10x",hist=False)
# sns.distplot( df_above10['at_count'] , color="green", label="AT% (avg)x",hist=False)
# sns.distplot( th_samples , color="blue", label="AT% Theoretical",hist=False)
# plt.title(tag_1+" AT distribution over sequences")
# plt.xlabel("Mean AT %")
# plt.ylabel("Proportion of reads")
# plt.legend()
# plt.savefig('ATrich_coverage_uniformity_test')
# plt.close()

# comparing with design file
#### just for design files
file1='/home/hhadmin/kits_comparison/v7_cre_comparision/in_CRE_corriell/DEMO.CREcoverage.mean.allx.fa.out'
df_1 = pd.read_csv(file1,sep='\t',header=None)
df_1_mod = getATCount(df_1)


sns.distplot( df_1_mod['at_count'] , color="green", label="cre_v2_design",hist=False)
sns.distplot( df_above10['at_count'] , color="blue", label=tag_1+"_average",hist=False)
sns.distplot( df_b10['at_count'] , color="red", label=tag_1+"_below_10x",hist=False)
plt.title("Comparison: AT distribution over sequences")
plt.xlabel("Mean AT %")
plt.ylabel("Proportion of reads")
plt.legend()
plt.savefig('Comparison_ATrich_uniformity_test')
plt.close()
