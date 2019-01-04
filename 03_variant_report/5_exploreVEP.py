import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib.pyplot import figure


def calcTiTvRatio(df):

    AtoG = df[ (df.REF=='A') & (df.ALT=='G') ].shape[0]
    GtoA = df[ (df.REF=='G') & (df.ALT=='A') ].shape[0]
    CtoT = df[ (df.REF=='C') & (df.ALT=='T') ].shape[0] ## very common
    TtoC = df[ (df.REF=='T') & (df.ALT=='C') ].shape[0]

    transition = AtoG + GtoA + CtoT + TtoC

    AtoC = df[ (df.REF=='A') & (df.ALT=='C') ].shape[0]
    CtoA = df[ (df.REF=='C') & (df.ALT=='A') ].shape[0]

    AtoT = df[ (df.REF=='A') & (df.ALT=='T') ].shape[0]
    TtoA = df[ (df.REF=='T') & (df.ALT=='A') ].shape[0]

    GtoC = df[ (df.REF=='G') & (df.ALT=='C') ].shape[0]
    CtoG = df[ (df.REF=='C') & (df.ALT=='G') ].shape[0]

    GtoT = df[ (df.REF=='G') & (df.ALT=='T') ].shape[0]
    TtoG = df[ (df.REF=='T') & (df.ALT=='G') ].shape[0]

    transversions = AtoC + CtoA + AtoT + TtoA + GtoC + CtoG + GtoT + TtoG

    return float(transition/transversions)


def calcIndelRatio(df):
    insertion = df[(df.REF.str.len()==1)& (df.ALT.str.len()>1)].shape[0]
    deletion = df[(df.REF.str.len()>1)& (df.ALT.str.len()==1)].shape[0]

    return float(insertion/deletion)





sample="COV_1_2_R"
df_vep = pd.read_csv("bcm_snp_indel_atlas.target.txt",sep='\t')
df_vep.columns=['CHROM','POS','REF','ALT']
print(calcTiTvRatio(df_vep))
print(calcIndelRatio(df_vep))
