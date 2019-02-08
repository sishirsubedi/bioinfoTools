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


def getFASTA(input,output):
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


DEPTH_CUTOFF=20
SAMPLE_COVERAGE_BED = sys.argv[1]
OUTPUT_DIR_SAMPLE=os.path.abspath(os.path.dirname(SAMPLE_COVERAGE_BED))
SAMPLE=os.path.basename(SAMPLE_COVERAGE_BED).split('.')[0]
print(SAMPLE)
print(OUTPUT_DIR_SAMPLE)


df_coverage_all = pd.read_csv(SAMPLE_COVERAGE_BED,header=None,sep='\t')
df_coverage_all.columns = ['chromosome','design_start','design_end','read_depth']
df_coverage_all['tag'] = SAMPLE

df_coverage_all.to_csv(OUTPUT_DIR_SAMPLE+"/"+SAMPLE+".CREcoverage.mean.allx.bed",sep='\t',header=False,index=False)

df_coverage_below = df_coverage_all[df_coverage_all['read_depth']<DEPTH_CUTOFF]
df_coverage_below.to_csv(OUTPUT_DIR_SAMPLE+"/"+SAMPLE+".CREcoverage.mean.b20x.bed",sep='\t',header=False,index=False)

df_coverage_above = df_coverage_all[(df_coverage_all['read_depth']>=DEPTH_CUTOFF)]
df_coverage_above.to_csv(OUTPUT_DIR_SAMPLE+"/"+SAMPLE+".CREcoverage.mean.above20x.bed",sep='\t',header=False,index=False)


###create fasta file for specific regions
getFASTA(OUTPUT_DIR_SAMPLE+"/"+SAMPLE+".CREcoverage.mean.allx.bed", OUTPUT_DIR_SAMPLE+"/"+SAMPLE+".CREcoverage.mean.allx.fa.out")
getFASTA(OUTPUT_DIR_SAMPLE+"/"+SAMPLE+".CREcoverage.mean.b20x.bed", OUTPUT_DIR_SAMPLE+"/"+SAMPLE+".CREcoverage.mean.b20x.fa.out")
getFASTA(OUTPUT_DIR_SAMPLE+"/"+SAMPLE+".CREcoverage.mean.above20x.bed", OUTPUT_DIR_SAMPLE+"/"+SAMPLE+".CREcoverage.mean.above20x.fa.out")

df_all = pd.read_csv( OUTPUT_DIR_SAMPLE+"/"+SAMPLE+".CREcoverage.mean.allx.fa.out",sep='\t',header=None)
getATCount(df_all)

df_below = pd.read_csv( OUTPUT_DIR_SAMPLE+"/"+SAMPLE+".CREcoverage.mean.b20x.fa.out",sep='\t',header=None)
getATCount(df_below)

df_above = pd.read_csv( OUTPUT_DIR_SAMPLE+"/"+SAMPLE+".CREcoverage.mean.above20x.fa.out",sep='\t',header=None)
getATCount(df_above)


plt.xlim(0,100)
sns.distplot( df_all['at_count'] , color="red", label="AT% (all)",hist=False)
sns.distplot( df_below['at_count'] , color="green", label="AT% (<20x)",hist=False)
sns.distplot( df_above['at_count'] , color="blue", label="AT% (>20x)",hist=False)
plt.title(SAMPLE+" AT distribution over sequences")
plt.xlabel("Mean AT %")
plt.ylabel("Proportion of reads")
plt.legend()
plt.savefig(OUTPUT_DIR_SAMPLE+"/"+SAMPLE+'_ATrich_coverage_uniformity_test')
plt.close()
