
import pandas as pd
import numpy as np
import argparse
import os
import re


def analyzeVCF(vcfs_directory,filter=False,out_file=None):

    vcfs = os.listdir(vcfs_directory)

    print("----VCF analysis ----")

    df_combine = pd.DataFrame()
    failed = 0
    for vcf_file in vcfs:
        try:
            df = pd.read_csv(vcfs_directory+vcf_file,sep="\t",skiprows=9)
            df["Strain"] = vcf_file.split(".")[0]
            df["Primary_Call"] = [x.split(";")[0].split("=")[1] for x in df.INFO.values]
            df = df[["Strain","POS","REF","ALT","QUAL","Primary_Call"]]
            df_combine = pd.concat([df_combine,df],axis=0)
        except:
            failed += 1
            print("failed number:"+str(failed)+"---"+vcf_file)

    print("Total strains failed i.e. could not read file : " + str(failed))
    print("Total strains with data : " + str(len(df_combine["Strain"].unique())))

    ## remove repeat variants
    df_combine["Strain"] = [x.replace("v3","").replace("combo","").replace("-r1","").replace("_r1","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-0","-") for x in df_combine["Strain"].values]

    df_combine.drop_duplicates(["Strain","POS","REF","ALT"],inplace=True)

    total_samples = len(df_combine["Strain"].unique())

    print("Total strains after removing repeats before vcf filter: " + str(total_samples))


    #### filters
    df_combine = df_combine[df_combine["ALT"] == df_combine["Primary_Call"]]

    #filter quality

    if filter:
        df_combine["QUAL"] = df_combine["QUAL"].astype(float)
        df_combine = df_combine[df_combine["QUAL"]>=50.0]

    if out_file != None:
        df_combine.to_csv(out_file,index=False)
    else:
        return df_combine

def getTotalNumberofVCFVariants(vcfs_directory,filter=False,out_file=None):

    vcfs = os.listdir(vcfs_directory)

    print("----VCF analysis ----")

    vcf_count =[]    
    failed = 0
    for vcf_file in vcfs:
        try:
            df = pd.read_csv(vcfs_directory+vcf_file,sep="\t",skiprows=9)
            vcf_count.append([vcf_file.split(".")[0],df.shape[0]])
        except:
            failed += 1
            print("failed number:"+str(failed)+"---"+vcf_file)

    ## remove repeat variants
    df_combine = pd.DataFrame(vcf_count)
    df_combine.columns = ["Strain","totalVariant"]
    df_combine["Strain"] = [x.replace("-r1","").replace("_r1","").replace("r1","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-r2","").replace("_r2","").replace("r2","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-r3","").replace("_r3","").replace("r3","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("v3","").replace("V3","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("_","-") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-0","-") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-R","") for x in df_combine["Strain"].values]

    df_combine.drop_duplicates(["Strain"],inplace=True)

    if out_file != None:
        df_combine.to_csv(out_file,index=False)
    else:
        return df_combine

def analyzePositionVariant(vcfs_directory,position,out_file):

    print ("processing.....vcfs")

    vcfs = os.listdir(vcfs_directory)
    combine = []
    failed = 0
    for vcf_file in vcfs:
        try:
            df = pd.read_csv(vcfs_directory+vcf_file,sep="\t",skiprows=9)
            df = df[['POS','REF','ALT','QUAL']]
            df.POS = df.POS.astype(int)

            #check for D614G mutation
            df = df[ ( (df.POS==position) & (df.REF=="A") & (df.ALT=="G") ) ]

            if df.shape[0] == 1:
                temp = []
                temp.append(vcf_file.split(".")[0])
                temp.append(df.iloc[0,0])
                temp.append(df.iloc[0,1])
                temp.append(df.iloc[0,2])
                temp.append(df.iloc[0,3])
                combine.append(temp)
            else:
                temp = []
                temp.append(vcf_file.split(".")[0])
                temp.append("NotFound")
                combine.append(temp)
        except:
            failed += 1
            print("failed number:"+str(failed)+"---"+vcf_file)
    df_combine = pd.DataFrame(combine,columns=["Strain","POS","REF","ALT","QUAL"])

    df_combine["Strain"] = [x.replace("-r1","").replace("_r1","").replace("r1","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-r2","").replace("_r2","").replace("r2","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-r3","").replace("_r3","").replace("r3","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("v3","").replace("V3","") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("_","-") for x in df_combine["Strain"].values]
    df_combine["Strain"] = [x.replace("-0","-") for x in df_combine["Strain"].values]

    print(df_combine.shape)
    ## if repeat then keep a strain with high quality variant
    df_combine.POS = df_combine.POS.astype(str)
    df_combine.sort_values(["POS","QUAL"],inplace=True)
    df_combine.drop_duplicates("Strain",inplace=True)
    print(df_combine.shape)

    if out_file != None:
        df_combine.to_csv(out_file,index=False)
    else:
        return df_combine

