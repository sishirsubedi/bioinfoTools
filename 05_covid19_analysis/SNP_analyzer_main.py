import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import seaborn as sns
from Bio import AlignIO, SeqIO, Entrez
import numpy as np
import os
import re
import subprocess

def bashCommunicator(command,output_expected=False):
    process = subprocess.Popen([command],shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print("Process failed %s %d %s %s" % (command,process.returncode, stdout, stderr))
    else:
        print("bash command completed.")
        if output_expected:
            return [x for x in stdout.split("\n")]

def analyzeVCF(vcfs_directory,region,region_name,out_dir,reported):

    vcfs = os.listdir(vcfs_directory)

    region_start = region[0]
    region_stop = region[1]

    df_combine = pd.DataFrame()
    failed = 0
    for vcf_file in vcfs:
        try:
            df = pd.read_csv(vcfs_directory+vcf_file,sep="\t",skiprows=9)
            df["STRAIN"] = vcf_file.split(".")[0]
            df["Primary_Call"] = [x.split(";")[0].split("=")[1] for x in df.INFO.values]
            df = df[["STRAIN","POS","REF","ALT","QUAL","Primary_Call"]]
            df_combine = pd.concat([df_combine,df],axis=0)
        except:
            failed += 1
            print("failed number:"+str(failed)+"---"+vcf_file)

    print("----VCF analysis ----")
    print("Total strains failed i.e. could not read file : " + str(failed))
    print("Total strains with data : " + str(len(df_combine["STRAIN"].unique())))

    ## remove repeat variants
    df_combine["STRAIN"] = [x.replace("v3","").replace("combo","").replace("-r1","").replace("_r1","") for x in df_combine["STRAIN"].values]

    df_combine.drop_duplicates(["STRAIN","POS","REF","ALT"],inplace=True)

    total_samples = len(df_combine["STRAIN"].unique())

    print("Total strains after removing repeats before vcf filter: " + str(total_samples))


    #### filters
    df_combine = df_combine[df_combine["ALT"] == df_combine["Primary_Call"]]

    #filter quality
    df_combine["QUAL"] = df_combine["QUAL"].astype(float)
    df_combine = df_combine[df_combine["QUAL"]>=50.0]


    df_combine = df_combine[( ( df_combine["POS"] >= region_start ) & ( df_combine["POS"] <= region_stop ) ) ]

    print("Total strains after removing repeats after vcf filter: " + str(len(df_combine["STRAIN"].unique())))


    ### remove N with >5%
    df_ncount = pd.read_excel(out_dir+"N_SNPcount.xlsx")
    df_ncount["Nproportion"] = df_ncount["Nproportion"].astype(float)
    df_ncount = df_ncount[df_ncount["Nproportion"]<=5.0]

    print("Total strains with <5% N from fasta data " + str(df_ncount.shape[0]))
    df_combine = df_combine[df_combine["STRAIN"].isin(df_ncount["Strain"].values)]
    print("Total strains after removing repeats after vcf filter: " + str(len(df_combine["STRAIN"].unique())))


    ref ={}
    ref[0]="ref"
    for indx,row in df_combine.iterrows():
        ref[row["POS"]]=row["REF"]

    df_ref = pd.DataFrame(
        data=ref,
        index=['REFERENCE'],
        columns=[x for x in ref.keys()])

    for strain in df_combine["STRAIN"].unique():
        for indx,row in df_combine[df_combine["STRAIN"]==strain].iterrows():
            df_ref.ix[row["STRAIN"],row["POS"]] = row["ALT"]

    mismatch =[]
    for column in df_ref:
        df_clean = df_ref[df_ref[column].notnull()]
        mismatch.append(sum(df_clean[column]!= df_clean.ix[0,column]))


    for column in df_ref:
        df_ref.ix[df_ref[column].isnull(),column] = df_ref.ix[0,column]


    df_ref = df_ref.T
    df_ref=df_ref.iloc[1:,:]
    df_ref["VCF_Mismatch_Count"] = mismatch[1:]

    sample_number = df_ref.shape[1]-1
    variant_present =[]
    variant_count =[]
    for indx,row in df_ref.iterrows():

        vc = str(row[0:sample_number].value_counts()).split("\n")
        variant_count.append(vc[:len(vc)-1])

        t = row[0:sample_number].unique()
        variant_present.append([x for x in t if x!= row[0]])


    df_ref.reset_index(inplace=True)
    variants_info =[]
    for indx,row in df_ref.iterrows():
        variants = {}
        for variant in row[1:sample_number].unique():
            strains = df_ref.ix[:,df_ref.loc[indx] == variant].columns
            variants[variant] = []
            variants[variant].append(len(strains))
            if "REFERENCE" in strains:
                variants[variant].append(["REFERENCE"])
            else:
                variants[variant].append(list(strains))
        variants_info.append(variants)

    df_ref["VCF_Variants"] = variant_present
    df_ref["VCF_Variants_Count"] = variant_count
    df_ref["VCF_Variants_Info"] = variants_info

    df_ref = df_ref[['index', 'REFERENCE','VCF_Mismatch_Count', 'VCF_Variants','VCF_Variants_Count','VCF_Variants_Info']]
    df_ref.columns = ['VCF_Genomic_Locus', 'VCF_Reference','VCF_Mismatch_Count', 'VCF_Variants','VCF_Variants_Count','VCF_Variants_Info']

    ### add clades column ###
    # df_clades = pd.read_csv("reference/clades.csv")
    # df_clades = df_clades[( ( df_clades["gene_position_start"] >= region_start ) & ( df_clades["gene_position_stop"] <= region_stop ) ) ]
    # df_ref["VCF_Clade"] =""
    # for indx,row in df_ref.iterrows():
    #     for indx2,row2 in df_clades.iterrows():
    #         if ( ( row["VCF_Genomic_Locus"] >= row2["gene_position_start"]) and ( row["VCF_Genomic_Locus"] <= row2["gene_position_stop"]) ):
    #             df_ref.ix[indx,"VCF_Clade"] = row2["clade"]


    ####add annotation from SNPEFF
    vcf_file = []
    chr_pos_ref_alt = []
    for indx,row in df_ref.iterrows():
        for indx2,variant in enumerate(row.VCF_Variants):

            vcf_file.append(["NC_045512.2",row.VCF_Genomic_Locus,".",row.VCF_Reference,variant,"100.0","PASS","INFO"])

            variant_occurence=0
            for variant_count in row.VCF_Variants_Count:
                if variant == variant_count[0]:
                    variant_occurence = int(variant_count.replace(variant,'').replace(' ',''))

            chr_pos_ref_alt.append([
            row.VCF_Genomic_Locus,row.VCF_Reference,row.VCF_Mismatch_Count,row.VCF_Variants,row.VCF_Variants_Count,row.VCF_Variants_Info,
            "NC_045512.2",row.VCF_Genomic_Locus,row.VCF_Reference,variant,variant_occurence])


    df_vcf_file = pd.DataFrame(vcf_file)
    df_vcf_file[1] = df_vcf_file[1].astype(int)
    df_vcf_file = df_vcf_file.drop_duplicates()
    df_vcf_file = df_vcf_file.sort_values(1)
    df_vcf_file.to_csv(out_dir+"covid.vcf",sep='\t',header=None,index=False)

    ###update df_ref in row wise split of variants
    new_columns =[x for x in df_ref.columns]
    new_columns.append("CHROM")
    new_columns.append("POS")
    new_columns.append("REF")
    new_columns.append("ALT")
    new_columns.append("ALT_count")
    df_ref2 = pd.DataFrame(chr_pos_ref_alt,columns=new_columns)

    cmd="rm -f %s" %(out_dir+"covid.vcf.snpeff")
    bashCommunicator(cmd)
    print("Delete old annotation file")

    cmd = "java -Xmx4g -jar /opt/snpeff/snpeff_covid/snpEff/snpEff.jar -v   NC_045512.2  %s > %s " %(out_dir+"covid.vcf",out_dir+"covid.vcf.snpeff" )
    bashCommunicator(cmd)
    print("Generating new annotation file")

    VEP_FILE=out_dir+"covid.vcf.snpeff"
    header=[]
    with open(VEP_FILE) as myfile:
        header = [next(myfile) for x in range(3)]
    columns_line=str(header[2].strip().split(':')[1]).split('|')

    df_snpeff = pd.read_csv(VEP_FILE,sep='\t',skiprows=5,header=None)
    df_snpeff.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    df_snpeff["HGMD_INFO"]= [x.split(";")[0] for x in df_snpeff.INFO]
    df_snpeff["VEP_INFO"]= [x.split(";")[1].split(',')[0] for x in df_snpeff.INFO]
    df_snpeff[columns_line] = df_snpeff['VEP_INFO'].str.split('|',expand=True)

    dfjoin = pd.merge(df_ref2,df_snpeff,how="left",on=["POS","REF","ALT"],indicator=True)


    gene_start = region_start+1
    if region_name == "S":
        gene_start = region_start-1

    dfjoin["Gene_Locus"] = [ int(x)-gene_start for x in dfjoin["VCF_Genomic_Locus"].values]
    dfjoin["Protein_AA_Locus"] = [ int(x)/3 for x in dfjoin["Gene_Locus"].values]
    dfjoin["Protein_AA_Locus"] = [ int(x) if x.is_integer() else int(x+1)  for x in dfjoin["Protein_AA_Locus"].values]

    df_reported = pd.read_csv(reported+region_name+"_curated.csv")
    curated = df_reported["Genomic locus"].values


    dfjoin["Status"] = ["Reported" if x in curated else "NotReported" for x in dfjoin["VCF_Genomic_Locus"].values]

    dfjoin.columns = [x.replace(" ","") for x in dfjoin.columns]

    dfjoin = dfjoin[dfjoin['Annotation']!="stop_gained"]
    ### confirm amino acid changes

    aa_dict={"Ala": "A",
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
    'Annotation', 'Annotation_Impact', 'HGVS.c', 'Gene_Locus','HGVS.p','HGVS.p.short' , 'Status',
    'VCF_Variants','VCF_Mismatch_Count', 'VCF_Variants_Count', 'VCF_Variants_Info']]


    dfjoin.columns = ['POS', 'REF', 'ALT','ALT_count',
                      'Annotation', 'Annotation_Impact', 'HGVS.c', 'Gene_Locus','HGVS.p','HGVS.p.short' , 'Status',
                      'ALL_Variants','Total_Mismatch_Count','Count_per_Variant', 'Variants_Info' ]

    dfjoin.sort_values(by=['POS'],inplace=True)
    dfjoin.to_excel(out_dir+"COVID_SNP_MAP_"+region_name+"_vcfs.xlsx",index=False)

def analyzeEntireGenome(alignment_file,out_dir):
    alignment = AlignIO.read(alignment_file,"fasta")
    filt = []
    for line in alignment:
        temp =[]
        temp.append(line.id)
        for nt in str(line.seq):
            temp.append(nt)
        filt.append(temp)

    df_genome = pd.DataFrame(filt)
    df_genome = df_genome[~df_genome[0].str.contains("EPI_ISL")]
    # df_genome = df_genome[(df_genome != 'n').all(axis=1)]
    df_genome = df_genome.drop(df_genome.columns[df_genome.iloc[0,:] =='-'],axis=1)
    df_genome.columns = range(df_genome.shape[1])

    total_samples = df_genome.shape[0]

    genome_mismatch =[]
    for column in df_genome:
        df_clean = df_genome[df_genome[column]!='n']
        genome_mismatch.append(sum(df_clean[column]!= df_clean.iloc[0,column]))

    pd.DataFrame(genome_mismatch).to_csv("BACKUP_GENOME_MISMATCH.csv")

    sns.lineplot(x=range(1,len(genome_mismatch),1), y=genome_mismatch[1:]) ##remove first sampleID to plot
    plt.title("MCoV samples (n="+total_samples+")")
    plt.ylabel("frequency of mismatch")
    plt.xlabel("genomic coordinates")
    plt.savefig(out_dir+"COVID_SNP_MAP_Genomic.png");
    plt.close()


    col = []
    col.append("id")
    for x in range(1,len(genome_mismatch),1): col.append(x)   ### +1 to match with ncbi indexing starting with 1
    df_genome.columns = col


    df_genome = df_genome.T
    df_genome.rename(columns=df_genome.iloc[0],inplace=True)
    df_genome=df_genome.iloc[1:,:]
    df_genome["mismatch_count"] = genome_mismatch[1:]

    df_genome = df_genome[df_genome["mismatch_count"]>=1]

    variants =[]
    for indx,row in df_genome.iterrows():
        t = row[0:df_genome.shape[1]-1].unique()
        variants.append([x for x in t if x!= row[0]])

    df_genome["variants"] = variants
    df_genome.reset_index(inplace=True)
    df_genome = df_genome[['index', 'MN908947','mismatch_count', 'variants']]

    df_annotation = pd.read_excel("nCOV_Variation_Annotation.xlsx")
    dfjoin = pd.merge(df_genome,df_annotation,how='left',left_on='index',right_on='Genome position')
    dfjoin.to_excel(out_dir+"COVID_SNP_MAP_Genomic_variants.xlsx",index=False)

def analyzeGenomicRegion(alignment_file,region,region_name,out_dir,reported):

    region_start = region[0]
    region_stop = region[1]

    filt = []
    #### read alignment files for reference
    alignment = AlignIO.read("data/4_15/Houston.4-14.clean.fa","fasta")
    for line in alignment:
        temp =[]
        temp.append(line.id)
        print(line.id)
        for nt in str(line.seq[region_start:region_stop]):
            temp.append(nt)
        filt.append(temp)
        break

    #### read alignment files
    alignment = AlignIO.read(alignment_file,"fasta")
    # filt = []
    for line in alignment:
        temp =[]
        temp.append(line.id)
        for nt in str(line.seq[region_start:region_stop]):
            temp.append(nt)
        filt.append(temp)

    df = pd.DataFrame(filt)
    # df[0]=[x.replace("V-",".").replace("-","").replace(".","-") for x in df[0].values]
    df.drop_duplicates(0,inplace=True)

    print("before removing ...")
    print(df.shape)

    ## remove 1255 and 1343 and two training runs
    df = df[df[0] != "MCoV-1255"]
    df = df[df[0] != "MCoV-1343"]
    df = df[df[0] != "MCoV-1483"]
    df = df[df[0] != "MCoV-743"]
    df = df[df[0] != "MCoV-1219"]

    print("after removing ...")
    print(df.shape)
    ### remove 320 strains from paper
    # df_320 = pd.read_csv("data/4_15/320_strains.csv")
    # match_strains =[]
    # for s in df_320["320_strains"].unique():match_strains.append(s)
    # df = df[~df[0].isin(match_strains)]
    print(df.shape)
    print(df.head())

    ###keep inhouse samples only
    df = df[~df[0].str.contains("EPI_ISL")]
    #### remove nt from samples where n is present
    # df = df[(df != 'n').all(axis=1)]
    ### remove any -
    df = df.drop(df.columns[df.iloc[0,:] =='-'],axis=1)


    mismatch =[]
    for column in df:
        if column ==0:
            mismatch.append(0)
        else:
            df_clean = df[df[column].isin(['a','t','g','c'])]
            mismatch.append(sum(df_clean[column]!= df_clean.iloc[0,column]))

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


    df2 = df.T
    df2.rename(columns=df2.iloc[0],inplace=True)
    df2=df2.iloc[1:,:]
    df2["mismatch_count"] = mismatch[1:]
    df2['gene']=region_name

    df2 = df2[df2["mismatch_count"]>0]

    variants =[]
    for indx,row in df2.iterrows():
        t = row[0:df2.shape[1]-1].unique()
        t = [str(x) for x in t]
        t2 = [x for x in t if  x!= row[0] ]
        t3 = [x for x in t2 if ('a' in x) or ('t' in x) or  ('g' in x) or ('c' in x)  ]
        variants.append(t3)

    df2["variants"] = variants

    sample_number = df2.shape[1]-1
    df2.reset_index(inplace=True)
    variants_info =[]
    for indx,row in df2.iterrows():
        variants = {}
        for variant in row[1:sample_number-1].unique():
            strains = df2.ix[:,df2.loc[indx] == variant].columns
            variants[variant] = []
            variants[variant].append(len(strains))
            if "MN908947" in strains:
                variants[variant].append(["REFERENCE"])
            elif "mismatch_count" in strains:
                continue
            else:
                variants[variant].append(list(strains))
        variants_info.append(variants)

    df2["variants_Info"] = variants_info

    df2 = df2[['index', 'MN908947','mismatch_count', 'gene', 'variants','variants_Info']]
    df2.columns = ['NTA_Genomic_Locus', 'NTA_Reference','NTA_Mismatch_Count', 'NTA_Gene', 'NTA_Variants','NTA_Variants_Info']

    ####add annotation

    ####add annotation from SNPEFF
    vcf_file = []
    chr_pos_ref_alt = []
    for indx,row in df2.iterrows():
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

    ###update df_ref in row wise split of variants
    new_columns =[x for x in df2.columns]
    new_columns.append("CHROM")
    new_columns.append("POS")
    new_columns.append("REF")
    new_columns.append("ALT")
    new_columns.append("ALT_count")
    df_ref = pd.DataFrame(chr_pos_ref_alt,columns=new_columns)


    cmd="rm -f %s" %(out_dir+"nta_covid.vcf.snpeff")
    bashCommunicator(cmd)
    print("Delete old annotation file")

    cmd = "java -Xmx4g -jar /opt/snpeff/snpeff_covid/snpEff/snpEff.jar -v   NC_045512.2  %s > %s " %(out_dir+"nta_covid.vcf",out_dir+"nta_covid.vcf.snpeff" )
    bashCommunicator(cmd)
    print("Generating new annotation file")

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

    dfjoin = pd.merge(df_ref,df_snpeff,how="left",on=["POS","REF","ALT"],indicator=True)

    gene_start = region_start+1
    if region_name == "S":
        gene_start = region_start-1

    dfjoin["Gene_Locus"] = [ int(x)-gene_start for x in dfjoin["NTA_Genomic_Locus"].values]
    dfjoin["Protein_AA_Locus"] = [ int(x)/3 for x in dfjoin["Gene_Locus"].values]
    dfjoin["Protein_AA_Locus"] = [ int(x) if x.is_integer() else int(x+1)  for x in dfjoin["Protein_AA_Locus"].values]

    df_reported = pd.read_csv(reported+region_name+"_curated.csv")
    curated = df_reported["Genomic locus"].values


    dfjoin["Status"] = ["Reported" if x in curated else "NotReported" for x in dfjoin["NTA_Genomic_Locus"].values]

    dfjoin.columns = [x.replace(" ","") for x in dfjoin.columns]

    dfjoin = dfjoin[dfjoin['Annotation']!="stop_gained"]
    dfjoin = dfjoin[dfjoin['Annotation']!="start_lost"]

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
    'Annotation', 'Annotation_Impact', 'HGVS.c', 'Gene_Locus','HGVS.p','HGVS.p.short' , 'Status',
    'NTA_Variants','NTA_Mismatch_Count', 'NTA_Variants_Info']]


    dfjoin.columns = ['POS', 'REF', 'ALT','ALT_count',
                      'Annotation', 'Annotation_Impact', 'HGVS.c', 'Gene_Locus','HGVS.p','HGVS.p.short' , 'Status',
                      'ALL_Variants','Total_Mismatch_Count','Variants_Info' ]

    dfjoin.sort_values(by=['POS'],inplace=True)
    dfjoin.to_excel(out_dir+"COVID_SNP_MAP_"+region_name+"_NTA.xlsx",index=False)

def analyzeProtein(alignment_file,protein_name,protein_reference,out_dir):

    filt = []

    ref = []
    ref.append(protein_name)
    for aa in protein_reference:
        ref.append(aa)

    filt.append(ref)

    alignment = AlignIO.read(alignment_file,"fasta")

    for line in alignment:
        temp =[]
        temp.append(line.id)
        for nt in str(line.seq):
            temp.append(nt)
        filt.append(temp)

    df = pd.DataFrame(filt)
    ###keep inhouse samples only
    df = df[~df[0].str.contains("EPI_ISL")]
    ### remove any -
    df = df.drop(df.columns[df.iloc[0,:] =='-'],axis=1)


    mismatch =[]
    for column in df:
        df_clean = df[df[column]!='-']
        df_clean = df_clean[df_clean[column]!='X']
        mismatch.append(sum(df_clean[column]!= df_clean.iloc[0,column]))


    sns.lineplot(x=range(1,len(mismatch),1), y=mismatch[1:]) ##remove first sampleID to plot
    plt.title("MCoV samples (n="+str(df.shape[0])+")")
    plt.ylabel("frequency of mismatch")
    plt.xlabel(protein_name)
    plt.savefig(out_dir+"COVID_SNP_MAP_"+protein_name+"_protein.png");
    plt.close()


    col = []
    col.append("id")
    for x in range(1,len(mismatch),1): col.append(x)   ### +1 to match with ncbi indexing starting with 1
    df.columns = col


    df2 = df.T
    df2.rename(columns=df2.iloc[0],inplace=True)
    df2=df2.iloc[1:,:]
    df2["mismatch_count"] = mismatch[1:]
    df2['gene']=protein_name

    df2 = df2[df2["mismatch_count"]>0]

    sample_number = df2.shape[1]-2
    variants =[]
    variant_count =[]
    for indx,row in df2.iterrows():

        vc = str(row[0:sample_number].value_counts()).split("\n")
        variant_count.append(vc[:len(vc)-1])

        t = row[0:sample_number].unique()
        variants.append([x for x in t if x!= row[0]])

    df2["variants"] = variants
    df2["variant_count"] = variant_count
    df2.reset_index(inplace=True)
    df2 = df2[['index', protein_name,'mismatch_count', 'gene', 'variants','variant_count']]
    df2.to_excel(out_dir+"COVID_SNP_MAP_"+protein_name+"_aa_variants.xlsx",index=False)

def analyzeFastaVariants(fasta_directory,vcfs_directory,out_dir):

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
    df_combine.to_excel(out_dir+"NCount.xlsx",index=False)

    print ("processing.....vcfs")

    vcfs = os.listdir(vcfs_directory)
    combine = []
    failed = 0
    for vcf_file in vcfs:
        try:
            df = pd.read_csv(vcfs_directory+vcf_file,sep="\t",skiprows=9)
            temp = []
            temp.append(vcf_file.split(".")[0])
            temp.append(df.shape[0])
            combine.append(temp)
        except:
            failed += 1
            print("failed number:"+str(failed)+"---"+vcf_file)
    df_combine2 = pd.DataFrame(combine,columns=["Strain","SNPcount"])
    df_combine2.to_excel(out_dir+"SNPcount.xlsx",index=False)

    dfjoin = pd.merge(df_combine, df_combine2,how="outer",on="Strain")
    dfjoin.to_excel(out_dir+"N_SNPcount.xlsx",index=False)
