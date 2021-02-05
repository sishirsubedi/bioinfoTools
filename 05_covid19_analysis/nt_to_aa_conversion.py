import pandas as pd
import argparse
import subprocess
import ast
from gen_utils.gen_io import bashCommunicator
import gen_utils.gen_covid as covid


def collect_mutations(df):
    print("Collecting mutations")
    mutations = []
    for indx,row in df.iterrows():
        for variant in row.mutation_info:
            if variant not in mutations:
                mutations.append(variant)
    return mutations

def mutations_to_vcf_file(mutations,out_variant):
    vcf_file = []
    for variant in mutations:
        vcf_file.append(["NC_045512.2",variant[1:len(variant)-1],".",variant[0].upper(),variant[len(variant)-1].upper(),"100.0","PASS","INFO"])

    df_vcf_file = pd.DataFrame(vcf_file)
    df_vcf_file[1] = df_vcf_file[1].astype(int)
    df_vcf_file = df_vcf_file.drop_duplicates()
    df_vcf_file = df_vcf_file.sort_values(1)
    df_vcf_file.to_csv(out_variant,sep='\t',header=None,index=False)

    cmd="rm -f %s" %(out_variant+".snpeff")
    bashCommunicator(cmd)
    print("Delete old annotation file")
    cmd = "java -Xmx4g -jar /opt/snpeff/snpeff_covid/snpEff/snpEff.jar -v   NC_045512.2  %s > %s " %(out_variant,out_variant+".snpeff" )
    bashCommunicator(cmd)
    print("Generating new annotation file")

def curationAfterAnnotation(VEP_FILE):
    print("Starting annotation curation")
    header=[]
    with open(VEP_FILE) as myfile:
        header = [next(myfile) for x in range(3)]
    columns_line=str(header[2].strip().split(':')[1]).replace(" ","").split('|')

    df_snpeff = pd.read_csv(VEP_FILE,sep='\t',skiprows=5,header=None)
    df_snpeff.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    df_snpeff["HGMD_INFO"]= [x.split(";")[0] for x in df_snpeff.INFO]
    df_snpeff["VEP_INFO"]= [x.split(";")[1].split(',')[0] for x in df_snpeff.INFO]
    df_snpeff[columns_line] = df_snpeff['VEP_INFO'].str.split('|',expand=True)

    ##get nt to aa dictionary
    df_snpeff.loc[(df_snpeff["HGVS.p"]==""),"HGVS.p"]="p.UTR"
    df_snpeff = df_snpeff[['POS', 'REF', 'ALT', 'Gene_Name', 'HGVS.p']]
    df_snpeff["mut-key"]=[x+str(y)+z for x,y,z in zip(df_snpeff.REF,df_snpeff.POS,df_snpeff.ALT)]
    df_snpeff["mut-value"]=[x+"-"+y.replace("p.","") for x,y in zip(df_snpeff.Gene_Name,df_snpeff["HGVS.p"])]
    return {x:y for x,y in zip(df_snpeff["mut-key"],df_snpeff["mut-value"])}

def add_aa_mutations(df,nt_aa_dict):
    
    aa_dict = covid.covid_basics.AA_TriToS

    print("adding aa mutations")
    strains = []
    for indx,row in df.iterrows():
        mutations =[]
        for variant in row.mutation_info:
                mutations.append(nt_aa_dict[variant])
        strains.append(mutations)
    return strains

def add_aa_mutations_column(df,wd):

    mutations = collect_mutations(df)
    mutations_to_vcf_file(mutations,wd+"2208_strains_genomewide_mutations.vcf")

    nt_aa_dict = curationAfterAnnotation(wd+"2208_strains_genomewide_mutations.vcf.snpeff")
    df["aa_info"] = add_aa_mutations(df,nt_aa_dict)
    return df



