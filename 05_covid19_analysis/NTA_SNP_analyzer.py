import pandas as pd
from Bio import AlignIO
import argparse
import subprocess
import ast

def bashCommunicator(command,output_expected=False):
    process = subprocess.Popen([command],shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print("Process failed %s %d %s %s" % (command,process.returncode, stdout, stderr))
    else:
        print("bash command completed.")
        if output_expected:
            return [x for x in stdout.split("\n")]

def sequenceFromAlignment(alignment_file,gene,coordinates,reference=False,adjustment=0):

    gene_start = int(coordinates.split(':')[0])
    gene_stop = int(coordinates.split(':')[1])

    sequences = []
    alignment = AlignIO.read(alignment_file,"fasta")

    for strain in alignment:
        strain_info =[]
        strain_info.append(strain.id)
        for nt in str(strain.seq[gene_start+adjustment:gene_stop+adjustment]):
            strain_info.append(nt.lower())
        sequences.append(strain_info)
        if reference:
            break

    return sequences

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

# def mismatchCount(df):
#     mismatch =[]
#     mismatch.append(0)
#     for column in df.columns[1:]:
#         df_clean = df[df[column].isin(['a','t','g','c'])]
#         mismatch.append(sum(df_clean[column]!= df_clean[df_clean["id"]=="MN908947"][column].values[0]))
#     return mismatch


def mismatchCount(df):
    mismatch =[]
    mismatch.append(0)
    for column in df.columns[1:]:
        ref_seq = df[df["id"]=="MN908947"][column].values[0]
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

def addVariantsInfo(df):
    sample_number = df.shape[1]-1
    df.reset_index(inplace=True)
    variants_info =[]
    for indx,row in df.iterrows():
        variants = {}
        for variant in row[1:sample_number-1].unique():
            strains = df.ix[:,df.loc[indx] == variant].columns
            variants[variant] = []
            variants[variant].append(len(strains))
            if "MN908947" in strains:
                variants[variant].append(["REFERENCE"])
            elif "mismatch_count" in strains:
                continue
            else:
                variants[variant].append(list(strains))
        variants_info.append(variants)
    return variants_info

def analyzeGenomicRegion(alignment_file,region,region_name):

    reference_alignment = "data/4_15/Houston.4-14.clean.fa"

    all_strains = sequenceFromAlignment(reference_alignment,region_name,region,reference=True)
    all_strains.extend(sequenceFromAlignment(alignment_file,region_name,region,adjustment=18))
    df = pd.DataFrame(all_strains)

    return df

def getVariantsTable(df,region_name):

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
    df2["variants_Info"] = addVariantsInfo(df2)

    df2 = df2[['index', 'MN908947','mismatch_count', 'gene', 'variants','variants_Info']]
    df2.columns = ['NTA_Genomic_Locus', 'NTA_Reference','NTA_Mismatch_Count', 'NTA_Gene', 'NTA_Variants','NTA_Variants_Info']

    return df2

def addAnnotation(df,out_dir):
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
    df_vcf_file.to_csv(out_dir+"nta_covid.vcf",sep='\t',header=None,index=False)

    cmd="rm -f %s" %(out_dir+"nta_covid.vcf.snpeff")
    bashCommunicator(cmd)
    print("Delete old annotation file")
    cmd = "java -Xmx4g -jar /opt/snpeff/snpeff_covid/snpEff/snpEff.jar -v   NC_045512.2  %s > %s " %(out_dir+"nta_covid.vcf",out_dir+"nta_covid.vcf.snpeff" )
    bashCommunicator(cmd)
    print("Generating new annotation file")

    ###update df_ref in row wise split of variants
    new_columns =[x for x in df.columns]
    new_columns.append("CHROM")
    new_columns.append("POS")
    new_columns.append("REF")
    new_columns.append("ALT")
    new_columns.append("ALT_count")
    df_ref = pd.DataFrame(chr_pos_ref_alt,columns=new_columns)
    return df_ref

def curationAfterAnnotation(df,region_name,out_dir):
    print("Starting annotation curation")

    VEP_FILE=out_dir+"nta_covid.vcf.snpeff"
    header=[]
    with open(VEP_FILE) as myfile:
        header = [next(myfile) for x in range(3)]
    columns_line=str(header[2].strip().split(':')[1]).split('|')

    df_snpeff = pd.read_csv(VEP_FILE,sep='\t',skiprows=5,header=None)
    df_snpeff.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    df_snpeff["HGMD_INFO"]= [x.split(";")[0] for x in df_snpeff.INFO]
    df_snpeff["VEP_INFO"]= [x.split(";")[1].split(',')[0] for x in df_snpeff.INFO]
    df_snpeff[columns_line] = df_snpeff['VEP_INFO'].str.split('|',expand=True)

    dfjoin = pd.merge(df,df_snpeff,how="left",on=["POS","REF","ALT"],indicator=True)

    region_start,region_stop=0,0
    if region_name == "nsp12":
        region_start = 13442
        region_stop = 16237
    elif region_name =="S":
        region_start = 21563
        region_stop = 25385

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

    aa_dict={
    "Ala": "A",
    "Asx": "B",
    "Cys": "C",
    "Asp": "D",
    "Glu": "E",
    "Phe": "F",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Xle": "J",
    "Lys": "K",
    "Leu": "L",
    "Met": "M",
    "Asn": "N",
    "Hyp": "O",
    "Pro": "P",
    "Gln": "Q",
    "Arg": "R",
    "Ser": "S",
    "Thr": "T",
    "Glp": "U",
    "Val": "V",
    "Trp": "W",
    "Ter": "X",
    "Tyr": "Y",
    "Glx": "Z"}

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

    dfjoin.to_excel(out_dir+"COVID_SNP_MAP_"+region_name+"_NTA.xlsx",index=False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="analyze SNPs from alignment file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment_file", help="alignment file")
    parser.add_argument("--region_name", help="protein name")
    parser.add_argument("--region", help="protein name")
    parser.add_argument("--out_dir", help="output file")
    args = parser.parse_args()

    print("processing--\
    alignment_file = "+args.alignment_file+"\
    region_name = "+args.region_name+"\
    region = "+args.region+"\
    out_dir= "+args.out_dir)

    alignment_file = args.alignment_file
    region_name = args.region_name
    region = args.region
    out_dir= args.out_dir

    df = analyzeGenomicRegion(alignment_file,region,region_name)
    df = filterMainAlignment(df)
    df = getVariantsTable(df,region_name)
    df_annot = addAnnotation(df,out_dir)
    curationAfterAnnotation(df_annot,region_name,out_dir)
