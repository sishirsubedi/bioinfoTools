# Based on this paper--
#
# SNP panel identification assay (SPIA): a
# genetic-based assay for the identification
# of cell lines
#
#  /opt/python3/bin/python3 /home/hhadmin/scripts/bioinfoTools/03_variant_report/9_compareSNPs.py /home/environments/ngs_test/exomeAnalysis/190531_NB552038_0012_AHLVT2BGX9/UNpaired/Exome8-N_S6/varscan/Exome8-N_S6.comb.vcf /home/environments/ngs_test/exomeAnalysis/190531_NB552038_0012_AHLVT2BGX9/UNpaired/Exome8-T_S7/varscan/Exome8-T_S7.comb.vcf
# RAW:
# 603770 783927 428533 608690 175237 29.023800453815195
# HQ:
# 72178 65698 13207 6727 58971 89.7607233096898
# normal tumor n-only t-only match match%


import sys
import pandas as pd
from scipy.stats import binom

def filterSNP(dfn):
    extra_cols=['GT','VF','DP','AD']
    dfn[extra_cols]=dfn['sample'].str.split(':',expand=True)

    depth_cols =['NDP','ADP']
    dfn[depth_cols]=dfn['AD'].str.split(',',expand=True)


    dfn[['VF','DP','NDP','ADP']]=dfn[['VF','DP','NDP','ADP']].apply(pd.to_numeric)
    dfn = dfn[['#CHROM', 'POS', 'REF','ALT','VF','DP','NDP','ADP']]

    ## only valid CHROM
    dfn = dfn[~dfn['#CHROM'].str.contains("_")]
    ##filter depth
    dfn = dfn[ ( dfn['DP']>9 ) & ( dfn['NDP']>9 ) & ( dfn['ADP']>9 )]
    ##filter frequency
    dfn = dfn[dfn['VF']>0.1]

    return dfn


def calcDistance(df1,df2):

    match = pd.merge(df1,df2,on=['#CHROM', 'POS', 'REF'],how='outer',indicator=True)

    vNSNPs = match[match['_merge']=='both'].shape[0]

    match = match[match['_merge']=='both']

    match = match[['#CHROM', 'POS', 'REF','ALT_x','ALT_y']]

    match['filter']=[1 if x[0]!=x[1] else 0 for x in zip(match['ALT_x'], match['ALT_y'])]

    print(binom.pmf(match[match['filter']==1].shape[0],vNSNPs,1/1000))

    return (match[match['filter']==1].shape[0] / vNSNPs)


def printRes(df1,df2):

    result=pd.merge(df1,df2,on=['#CHROM', 'POS', 'REF', 'ALT'],how='outer',indicator=True)
    normal = result[result['_merge']=='left_only'].shape[0]
    tumor = result[result['_merge']=='right_only'].shape[0]
    match = result[result['_merge']=='both'].shape[0]

    match_p=0
    if df1.shape[0] > df2.shape[0] :
        match_p = match/df2.shape[0] * 100
    else :
        match_p = match/df1.shape[0] * 100

    print (df1.shape[0],',', df2.shape[0],',', normal,',', tumor,',', match,',', match_p)

    return (normal+tumor, match)



snpf1=sys.argv[1]
snpf2=sys.argv[2]

df1=pd.read_csv(snpf1,sep='\t',skiprows=1)
df2=pd.read_csv(snpf2,sep='\t',skiprows=1)

print("RAW:")
printRes(df1,df2)

df1_f1 = filterSNP(df1)
df2_f1 = filterSNP(df2)

print("HQ:")
df_match = printRes(df1_f1,df2_f1)

print("Distance:",calcDistance(df1_f1,df2_f1))
