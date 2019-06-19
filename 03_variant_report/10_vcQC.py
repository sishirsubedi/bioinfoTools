import sys
import pandas as pd

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

vc1=sys.argv[1]
df1=pd.read_csv(vc1,sep='\t')
print("titv ratio:",calcTiTvRatio(df1))
# print("indel ratio:",calcIndelRatio(df1))

# example:
#  /opt/python3/bin/python3
#        /home/hhadmin/scripts/bioinfoTools/03_variant_report/10_vcQC.py
#       /home/environments/ngs_test/exomeAnalysis/190612_NB552038_0014_AHLTKKBGX9/Paired/Exome22-T_S5_Exome22-N_S4/Exome22-T_S5_Exome22-N_S4.variantcallers.combinev2.10_10_2.vep.parse.txt
