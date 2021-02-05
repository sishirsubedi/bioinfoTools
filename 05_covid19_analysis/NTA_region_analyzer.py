import pandas as pd
from Bio import AlignIO
import argparse
import ast
import collections, numpy
from alignment.alignment_tools import analyzeGenomicRegionwithReference

def filterMainAlignment(df):

    print("before removing repeat markers..."+str(df.shape))

    ### remove repeat markers
    df[0] = [x.replace("-0","-") for x in df[0].values]
    df[0] = [x.replace("-R","") for x in df[0].values]
    df[0] = [x.replace("-r1","").replace("_r1","").replace("r1","") for x in df[0].values]
    df[0] = [x.replace("-r2","").replace("_r2","").replace("r2","") for x in df[0].values]
    df[0] = [x.replace("-r3","").replace("_r3","").replace("r3","") for x in df[0].values]
    df[0] = [x.replace("v3","").replace("V3","") for x in df[0].values]
    df[0] = [x.replace("_","-") for x in df[0].values]

    print("before removing duplicates..."+str(df.shape))
    df.drop_duplicates(0,inplace=True)

    print("before removing training runs..."+str(df.shape))
    ## remove training runs
    ### more noise data from 5x depth
    remove_index = df[df[0].isin(["MCoV-1255","MCoV-1343","MCoV-743",
    "MCoV-1219","MCoV-12783","MCoV-12765","MCoV-12805","MCoV-12123",
    "MCoV-12291","MCoV-12788","MCoV-12824","MCoV-12867","MCoV-11867"])].index
    df.drop(remove_index,inplace=True)

    print("before removing epi and nonreference columns..."+str(df.shape))
    ###keep inhouse samples only
    df = df[~df[0].str.contains("EPI_ISL")]
    df = df.drop(df.columns[df.iloc[0,:] =='-'],axis=1)

    print("after completing filter...")
    print(df.shape)
    print(df.head())

    return df


def updateCoordinates(df,region_name):

    if region_name == "nsp12":
        region_start = 13442
        region_stop = 16237
    elif region_name =="S":
        region_start = 21563
        region_stop = 25385

    col = []
    col.append("id")
    for x in range(region_start,region_stop,1): col.append(x)   ### +1 to match with ncbi indexing starting with 1
    df.columns = col

    return df



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="analyze SNPs from alignment file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment_file", help="alignment file")
    parser.add_argument("--protein", help="protein name")
    parser.add_argument("--out_file", help="output file")
    args = parser.parse_args()

    alignment_file = args.alignment_file
    protein = args.protein
    out_file= args.out_file

    print(alignment_file)
    print(protein)
    print(out_file)


    df = analyzeGenomicRegionwithReference(alignment_file,protein)
    df = filterMainAlignment(df)
    df = updateCoordinates(df,protein)

    df.to_csv(out_file,index=False)