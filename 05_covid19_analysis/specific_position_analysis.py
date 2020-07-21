""" analyze specific variant using vcf, consensus fasta, and nuleotide alignment

##### example:

 /opt/python3/bin/python3 checkD614GNpercentStrains.py \
 --request "get-strains" \
 --alignment "/home/tmhsxs240/COVID_19/data/6_24/Houston.July1.clean--RedundantMRN.fa" \
 --output "/home/tmhsxs240/COVID_19/data/Ns_Strain_Comparison/D614G_Npercent_analysis.xlsx"

"""

import pandas as pd
import numpy as np
import argparse
import os
import re
from Bio import AlignIO, SeqIO, Entrez
import alignment_analysis
import fasta_analysis
import vcf_analysis


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign mutation status, sequence quality, and alignment status to strain",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--position", help="position for mutation of interest")
    parser.add_argument("--vcf", help="vcf directory")
    parser.add_argument("--fasta", help="fasta directory")
    parser.add_argument("--alignment", help="alignment file")
    parser.add_argument("--output", help="output file")
    args = parser.parse_args()


    df_vcfs = vcf_analysis.analyzePositionVariant(args.vcfs,int(args.position))
    df_npercent = fasta_analysis.getSequenceQuality(args.consensus)
    df_alignment = alignment_analysis.getStrains(args.alignment)

    dfjoin = pd.merge(df_vcfs,df_npercent,how="left",on="Strain")
    dfjoin["Present in NTA (1484 strains)"] = [1 if x in df_alignment.Strain.values else 0 for x in dfjoin.Strain.values]
    dfjoin.to_excel(args.output,index=False)
