import pandas as pd
import argparse
import subprocess
import ast
import gen_utils.gen_covid as covid

def curationAfterAnnotation(variant_base,VEP_FILE,region_name,out_file):
    print("Starting annotation curation")
    header=[]
    with open(VEP_FILE) as myfile:
        header = [next(myfile) for x in range(3)]
    columns_line=str(header[2].strip().split(':')[1]).split('|')

    df_snpeff = pd.read_csv(VEP_FILE,sep='\t',skiprows=5,header=None)
    df_snpeff.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    df_snpeff["HGMD_INFO"]= [x.split(";")[0] for x in df_snpeff.INFO]
    df_snpeff["VEP_INFO"]= [x.split(";")[1].split(',')[0] for x in df_snpeff.INFO]
    df_snpeff[columns_line] = df_snpeff['VEP_INFO'].str.split('|',expand=True)

    df = pd.read_csv(variant_base)
    dfjoin = pd.merge(df,df_snpeff,how="left",on=["POS","REF","ALT"],indicator=True)

    region_start,region_stop=0,0
    if region_name == "nsp12":
        region_start = covid.covid_basics.nsp12_region_start
        region_stop = covid.covid_basics.nsp12_region_stop
    elif region_name =="S":
        region_start = covid.covid_basics.s_region_start
        region_stop = covid.covid_basics.s_region_stop

    gene_start = region_start+1
    if region_name == "S":
        gene_start = region_start-1

    dfjoin["Gene_Locus"] = [ int(x)-gene_start for x in dfjoin["NTA_Genomic_Locus"].values]
    dfjoin["Protein_AA_Locus"] = [ int(x)/3 for x in dfjoin["Gene_Locus"].values]
    dfjoin["Protein_AA_Locus"] = [ int(x) if x.is_integer() else int(x+1)  for x in dfjoin["Protein_AA_Locus"].values]

    dfjoin.columns = [x.replace(" ","") for x in dfjoin.columns]

    dfjoin = dfjoin[dfjoin['Annotation']!="stop_gained"]
    dfjoin = dfjoin[dfjoin['Annotation']!="start_lost"]
    dfjoin = dfjoin[dfjoin['Annotation']!="initiator_codon_variant"]
    dfjoin = dfjoin[dfjoin['Annotation']!="stop_lost&splice_region_variant"]

    ### confirm amino acid changes

    aa_dict= covid.covid_basics.AA_TriToS

    print(dfjoin.head())
    dfjoin["aa_auto1"]=[ x.split(".")[1] for x in dfjoin["HGVS.p"]]
    dfjoin["aa_auto2"]=[str(x[0:3])+'-'+str(x[len(x)-3:]) for x in dfjoin["aa_auto1"]]
    dfjoin["aa_auto3"]=[str(x[3:len(x)-3:]) for x in dfjoin["aa_auto1"]]

    if region_name == 'S':
        dfjoin["HGVS.p.short"]=[ aa_dict[x.split("-")[0]]+y+aa_dict[x.split("-")[1]] for x,y in zip(dfjoin["aa_auto2"],dfjoin["aa_auto3"])]
    elif region_name == 'nsp12':
        dfjoin["HGVS.p.short"]=[ aa_dict[x.split("-")[0]]+str(y)+aa_dict[x.split("-")[1]] for x,y in zip(dfjoin["aa_auto2"],dfjoin["Protein_AA_Locus"])]

    dfjoin["Gene_Locus"]=[ x.split(">")[0][len(x.split(">")[0])-1]+str(y)+ x.split(">")[1] for x,y in zip(dfjoin["HGVS.c"],dfjoin["Gene_Locus"])]


    dfjoin = dfjoin[['POS', 'REF', 'ALT','ALT_count',
    'Annotation', 'Annotation_Impact', 'HGVS.c', 'Gene_Locus','HGVS.p','HGVS.p.short' ,
    'NTA_Variants','NTA_Mismatch_Count', 'NTA_Variants_Info']]


    dfjoin.columns = ['POS', 'REF', 'ALT','ALT_count',
                      'Annotation', 'Annotation_Impact', 'HGVS.c', 'Gene_Locus','HGVS.p','HGVS.p.short' ,
                      'ALL_Variants','Total_Mismatch_Count','Variants_Info' ]

    dfjoin.sort_values(by=['POS'],inplace=True)

    dfjoin.to_excel(out_file,index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="analyze SNPs from alignment file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--variant_base", help="alignment file")
    parser.add_argument("--vep_file", help="vep out file")
    parser.add_argument("--region", help="protein name")
    parser.add_argument("--out_file", help="output file")
    args = parser.parse_args()

    curationAfterAnnotation(args.variant_base,args.vep_file,args.region,args.out_file)
