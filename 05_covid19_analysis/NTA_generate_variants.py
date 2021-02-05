import pandas as pd
from Bio import AlignIO
import argparse
import subprocess
import ast
import gen_utils.gen_covid as covid

def mismatchCount(df):
    mismatch =[]
    mismatch.append(0)
    for column in df.columns[1:]:
        ref_seq = df[df["id"]==covid.covid_basics.reference][column].values[0]
        var_table = df[column].value_counts()
        var_sum=0
        for v in var_table.keys():
            if v!=ref_seq and v in ["a","t","g","c"]:
                var_sum += var_table[v]
        mismatch.append(var_sum)
    return mismatch

def addVariants(df):
    variants =[]
    for indx,row in df.iterrows():
        t = row[0:df.shape[1]-1].unique()
        t = [str(x) for x in t]
        t2 = [x for x in t if  x!= row[0] ]
        t3 = [x for x in t2 if ('a' in x) or ('t' in x) or  ('g' in x) or ('c' in x)  ]
        variants.append(t3)
    return variants

def addVariantsInfo(df,mode):
    sample_number = df.shape[1]-1
    df.reset_index(inplace=True)
    variants_info =[]
    for indx,row in df.iterrows():
        variants = {}
        for variant in row[1:sample_number-1].unique():
            strains = df.loc[:,df.loc[indx] == variant].columns
            variants[variant] = []
            variants[variant].append(len(strains))
            if covid.covid_basics.reference in strains:
                variants[variant].append(["REFERENCE"])
            elif "mismatch_count" in strains:
                continue
            else:
                if mode == "all":
                    variants[variant].append(list(strains))
                elif mode == "single":
                    variants[variant].append(list(strains)[0])
        variants_info.append(variants)
    return variants_info

def addVariantsInfo_v2(data_mat):
    import collections
    variants_info =[]
    for indx in range(len(data_mat)):
        variants_info.append(dict(collections.Counter(data_mat[indx])))
    return variants_info

def getVariantsTable(alignment_file,region_name,strain_info_mode):

    df = pd.read_csv(alignment_file)

    print("Calculating mismatch")

    mismatch = mismatchCount(df)

    print("Completed calculating mismatch")

    df2 = df.T
    df2.rename(columns=df2.iloc[0],inplace=True)
    df2=df2.iloc[1:,:]
    df2["mismatch_count"] = mismatch[1:]
    df2['gene']=region_name
    df2 = df2[df2["mismatch_count"]>0]

    print("Calculating variants")
    df2["variants"] = addVariants(df2)

    print("Calculating strains with variants")
    df2["variants_Info"] = addVariantsInfo(df2,strain_info_mode)


    df2 = df2[['index', 'MN908947','mismatch_count', 'gene', 'variants','variants_Info']]
    df2.columns = ['NTA_Genomic_Locus', 'NTA_Reference','NTA_Mismatch_Count', 'NTA_Gene', 'NTA_Variants','NTA_Variants_Info']

    return df2

def addAnnotation(df,out_variant,out_variant_base):
    print("Starting annotation")
    vcf_file = []
    chr_pos_ref_alt = []
    for indx,row in df.iterrows():
        for indx2,variant in enumerate(row.NTA_Variants):

            vcf_file.append(["NC_045512.2",row.NTA_Genomic_Locus,".",row.NTA_Reference.upper(),variant.upper(),"100.0","PASS","INFO"])

            variant_occurence=0
            for key,data in row.NTA_Variants_Info.items():
                if variant == key:
                    variant_occurence = int(data[0])

            chr_pos_ref_alt.append([
            row.NTA_Genomic_Locus,row.NTA_Reference.upper(),row.NTA_Mismatch_Count,row.NTA_Gene,row.NTA_Variants,row.NTA_Variants_Info,
            "NC_045512.2",row.NTA_Genomic_Locus,row.NTA_Reference.upper(),variant.upper(),variant_occurence])


    df_vcf_file = pd.DataFrame(vcf_file)
    df_vcf_file[1] = df_vcf_file[1].astype(int)
    df_vcf_file = df_vcf_file.drop_duplicates()
    df_vcf_file = df_vcf_file.sort_values(1)
    df_vcf_file.to_csv(out_variant,sep='\t',header=None,index=False)

    ###update df_ref in row wise split of variants
    new_columns =[x for x in df.columns]
    new_columns.append("CHROM")
    new_columns.append("POS")
    new_columns.append("REF")
    new_columns.append("ALT")
    new_columns.append("ALT_count")
    df_ref = pd.DataFrame(chr_pos_ref_alt,columns=new_columns)
    df_ref.to_csv(out_variant_base,index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="analyze SNPs from alignment file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment_file", help="alignment file")
    parser.add_argument("--protein", help="protein name")
    parser.add_argument("--mode", help="variants mode name")
    parser.add_argument("--out_variant", help="output file")
    parser.add_argument("--out_variant_base", help="output file")
    args = parser.parse_args()

    alignment_file = args.alignment_file
    protein = args.protein
    strain_info_mode = args.mode
    out_variant= args.out_variant
    out_variant_base= args.out_variant_base

    print(alignment_file)
    print(protein)
    print(strain_info_mode)
    print(out_variant)
    print(out_variant_base)


    df2 = getVariantsTable(alignment_file,protein,strain_info_mode)
    addAnnotation(df2,out_variant,out_variant_base)
