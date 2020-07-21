""" analyze vcf files

##### example:

 /opt/python3/bin/python3 vcf_analysis.py \
 --request "get-strains" \
 --alignment "/home/tmhsxs240/COVID_19/data/6_24/6_24_vcfs/" \
 --output "/home/tmhsxs240/COVID_19/data/Ns_Strain_Comparison/vcfs_analysis.csv"

"""

import pandas as pd
import numpy as np
import argparse
import os
import re
from Bio import AlignIO, SeqIO, Entrez

def getSequenceQuality(fasta_directory,out_file=None):

    print ("processing.....fastas")
    fastas = os.listdir(fasta_directory)

    combine = []
    failed = 0
    for fasta_file in fastas:
        try:
            df = pd.read_csv(fasta_directory+fasta_file,skiprows=1,header=None)
            temp = []
            temp.append(fasta_file.split(".")[0])
            sequence = str(df.loc[0:0].values[0])
            ncount = len(re.findall("N",sequence))
            temp.append(len(sequence))
            temp.append(ncount)
            temp.append( (ncount/len(sequence) *100 ))
            combine.append(temp)

        except:
            failed += 1
            print("failed number:"+str(failed)+"---"+fasta_file.split(".")[0])

    df_combine = pd.DataFrame(combine,columns=["Strain","Sequence","Ncount","Nproportion"])


    df_combine["Strain"] = [x.replace("-r1","").replace("_r1","").replace("r1","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-r2","").replace("_r2","").replace("r2","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-r3","").replace("_r3","").replace("r3","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("v3","").replace("V3","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("_","-") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-0","-") for x in df_combine["Strain"].values]

    df_combine.sort_values(["Nproportion"],inplace=True)
    df_combine.drop_duplicates("Strain",inplace=True)

    if out_file != None:
        df_combine.to_csv(out_file,index=False)
    else:
        return df_combine

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign mutation status, sequence quality, and alignment status to strain",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--request", help="task to perform")
    parser.add_argument("--vcf", help="vcfs file")
    parser.add_argument("--output", help="output file")
    args = parser.parse_args()

    if args.request == "get-strains":
        getStrains(args.alignment,args.output)
    elif args.request == "mutations-table":
        tableGenomicMutations(args.alignment,args.output)
