import pandas as pd
import argparse
import subprocess
import ast
from Bio import SeqIO
from gen_utils.gen_io import bashCommunicator
import gen_utils.gen_covid as covid


def calculateMutations(var_list):
    mut_counter = 0
    mut_list = []
    for var in var_list:
        if var[1] != var[2]:
            mut_counter += 1
            mut_list.append(var[1].upper()+str(var[0])+var[2].upper())
    return (mut_counter,mut_list)

def genomeWideMutationTable(df):
    REFERENCE = "MN908947"
    df2 = df.T
    df2.rename(columns=df2.iloc[0],inplace=True)
    df2=df2.iloc[1:,:]
    mismatch =[]
    for ci,column in enumerate(df2.columns[1:]):
        df_sel = df2[[REFERENCE,column]][df2[column].isin(['a','t','g','c','-'])]
        mut_counter,mut_list = calculateMutations(zip(df_sel.index,df_sel[REFERENCE],df_sel[column]))
        mismatch.append([column,mut_counter,mut_list])
    df_mutationTable = pd.DataFrame(mismatch)
    df_mutationTable.columns = ["strain","total_mutations","nucleotide_change"]
    return df_mutationTable


def collect_deletions(df):
    print("Collecting mutations")
    deletions = []
    for indx,row in df.iterrows():
        del_mutations = []
        for variant in row["nucleotide_change"]:
            if '-' in variant:##for deletion
                if len(variant)>5 and int(variant[1:len(variant)-1])<29675:##end of CDS
                    del_mutations.append(variant)
        deletions.append([row['strain'],del_mutations])
    return deletions

def deletions_for_vcf(deletions_strainwise):
    deletions_dict_strainwise = []
    for del_entry in deletions_strainwise:
        strain=del_entry[0]
        deletions = del_entry[1]
        deletions.append("last")
        del_dict = {}
        indx = 0
        d = deletions[indx]
        if d != "last":
            current_coordinate = int(d[1:len(d)-1])
            current_key= current_coordinate
            del_dict[current_coordinate] = str(d[0])
            while deletions[indx+1] != "last":
                next_d = deletions[indx+1]
                next_coordinate = int(next_d[1:len(next_d)-1])
                if  next_coordinate == current_coordinate +1:
                    d = deletions[indx+1]
                    del_dict[current_key] += str(d[0])
                    indx = indx+1
                    current_coordinate = current_coordinate +1
                else:
                    del_dict[next_coordinate] = str(next_d[0])
                    current_key = next_coordinate
                    current_coordinate = next_coordinate
                    indx = indx+1
        deletions_dict_strainwise.append([strain,del_dict])
    return deletions_dict_strainwise


def add_deletions_to_vcf(df,out_variant):

    del_list = collect_deletions(df)

    del_list_withstrain = deletions_for_vcf(del_list)

    reference_file = "/storage/scratch/covid/container/reference/MN908947.fasta"
    reference = SeqIO.read(reference_file,"fasta")

    vcf_file = []
    
    for del_entry in del_list_withstrain:
        strain = del_entry[0]
        del_dict = del_entry[1]
        for coordinate in del_dict:
            deleted_nt = del_dict[coordinate]
            position = coordinate -1
            ref = reference[position-1] +deleted_nt
            alt = reference[position-1]
            vcf_file.append(["NC_045512.2",position,".",ref,alt,"100.0",strain,"INFO"])

    df_vcf_file = pd.DataFrame(vcf_file)
    df_vcf_file[1] = df_vcf_file[1].astype(int)
    df_vcf_file = df_vcf_file.drop_duplicates()
    df_vcf_file = df_vcf_file.sort_values(1)
    df_vcf_file.to_csv(out_variant,sep='\t',header=None,index=False)

    cmd="rm -f %s" %(out_variant+".snpeff")
    bashCommunicator(cmd)
    print("Delete old Deletion annotation file")
    cmd = "java -Xmx4g -jar /storage/apps/opt/snpEff/snpEff_covid/snpEff.jar -v   NC_045512.2  %s > %s " %(out_variant,out_variant+".snpeff" )
    bashCommunicator(cmd)
    print("Generating new Deletion annotation file")

def collect_mutations(df):
    print("Collecting mutations")
    mutations = []
    for indx,row in df.iterrows():
        for variant in row["nucleotide_change"]:
            if '-' in variant:##for deletion
                continue
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
    cmd = "java -Xmx4g -jar /storage/apps/opt/snpEff/snpEff_covid/snpEff.jar -v   NC_045512.2  %s > %s " %(out_variant,out_variant+".snpeff" )
    bashCommunicator(cmd)
    print("Generating new annotation file")


def readVEP(VEP_FILE):
    header=[]
    with open(VEP_FILE) as myfile:
        header = [next(myfile) for x in range(3)]
    columns_line=str(header[2].strip().split(':')[1]).replace(" ","").split('|')

    df_snpeff = pd.read_csv(VEP_FILE,sep='\t',skiprows=5,header=None)
    df_snpeff.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    df_snpeff["HGMD_INFO"]= [x.split(";")[0] for x in df_snpeff.INFO]
    df_snpeff["VEP_INFO"]= [x.split(";")[1].split(',')[0] for x in df_snpeff.INFO]
    df_snpeff[columns_line] = df_snpeff['VEP_INFO'].str.split('|',expand=True)

    return df_snpeff

def curationAfterAnnotation(VEP_FILE):
    print("Starting annotation curation")
    df_snpeff = readVEP(VEP_FILE)
    ##get nt to aa dictionary
    ##b/c of issue with UTR and empty annotation mark mutations with no protein change
    df_snpeff.loc[(df_snpeff["HGVS.p"]==""),"HGVS.p"]="p.UNKNOWN"

    df_snpeff = df_snpeff[['POS', 'REF', 'ALT', 'Gene_Name', 'HGVS.p']]
    df_snpeff["mut-key"]=[x+str(y)+z for x,y,z in zip(df_snpeff.REF,df_snpeff.POS,df_snpeff.ALT)]
    df_snpeff["mut-value"]=[x+"-"+y.replace("p.","") for x,y in zip(df_snpeff.Gene_Name,df_snpeff["HGVS.p"])]
    return {x:y for x,y in zip(df_snpeff["mut-key"],df_snpeff["mut-value"])}

def add_mutation_pair(df,nt_aa_dict):

    aa_dict = covid.covid_basics.AA_TriToS

    print("adding nt-aa pair mutations")
    strains = []
    for indx,row in df.iterrows():
        mutations =[]
        for variant in row["nucleotide_change"]:

            if '-' in variant:
                continue

            mutation = nt_aa_dict[variant]

            if 'UNKNOWN' in mutation:
                continue

            aa_change = mutation.split("-")[1]
            ref = aa_change[0:3]
            position = "".join([x for x in aa_change[3:] if x.isdigit()])

            alt = aa_change.replace(ref,"").replace(position,"")
            if alt =="": alt=ref ## no aa change for syn

            if alt in aa_dict.keys():
                alt= aa_dict[alt]
            mutations.append(variant+":"+mutation.split("-")[0]+"-"+aa_dict[ref]+position+alt)
        strains.append(mutations)
    return strains

def add_deletion_column(df,VEP_FILE):
    print("Starting Deletion annotation curation")
    df_snpeff = readVEP(VEP_FILE)
    df_snpeff =  df_snpeff[df_snpeff["HGVS.p"]!=""]

    aa_dict = covid.covid_basics.AA_TriToS
    del_col = []
    for indx,row in df.iterrows():
        dels = df_snpeff[df_snpeff["FILTER"]==row['strain']][["Gene_Name","HGVS.p"]].values

        if len(dels) == 0:
            del_col.append("")
        else:
            dnew=[]
            for d in dels:
                
                gene = d[0]
                aa_change = d[1]

                if "p." in aa_change:
                    aa_change = aa_change.replace("p.","")
                for aa in aa_dict:
                    if aa in aa_change:
                        aa_change = aa_change.replace(aa,aa_dict[aa])
                dnew.append(gene+":"+aa_change)
            del_col.append(dnew)

    return del_col


def add_aa_mutations_column(df,wd,tag):

    mutations = collect_mutations(df)
    mutations_to_vcf_file(mutations,wd+"strains_genomewide_mutations_"+tag+".vcf")

    nt_aa_dict = curationAfterAnnotation(wd+"strains_genomewide_mutations_"+tag+".vcf.snpeff")
    df["nt_aa_pair"] = add_mutation_pair(df,nt_aa_dict)

    return df

def add_aa_deletions_column(df,wd,tag):
    add_deletions_to_vcf(df,wd+"strains_genomewide_mutations_"+tag+".del.vcf")
    df["deletions"] = add_deletion_column(df,wd+"strains_genomewide_mutations_"+tag+".del.vcf.snpeff")
    return df

def generate_snp_table(df_muttable,db,out_dir,tag):

    df_db = pd.read_csv(db)

    sel_runs = sorted([x for x in df_db.run_id_seq.unique() if "low_quality" not in x])

    print(sel_runs)
    print(df_muttable.head(3))
    per_run_dfs = []

    for runid in sel_runs:

        print("current run..."+runid)

        sel_strains = list(df_db[df_db.run_id_seq==runid]["MCoVNumber"])
        df_run = df_muttable[df_muttable.strain.isin(sel_strains)]
        df_run = df_run[["strain","nt_aa_pair"]]

        if df_run.shape[0] == 0:
            continue

        spike_mutations = []
        for indx,row in df_run.iterrows():
            for mutation in row.nt_aa_pair:
                if mutation.split(':')[1].split('-')[0] =="S":

                    nt_change = mutation.split(':')[0]
                    gene = mutation.split(':')[1].split('-')[0]
                    aa_change = mutation.split(':')[1].split('-')[1]

                    if "*" in aa_change or "?" in aa_change:
                        continue

                    nt_ref = nt_change[0]
                    nt_alt = nt_change[len(nt_change)-1]
                    nt_position= nt_change.replace(nt_ref,"").replace(nt_alt,"")

                    aa_ref = aa_change[0]
                    aa_alt = aa_change[len(aa_change)-1]
                    aa_position= int("".join([x for x in aa_change if x.isdigit()]))

                    if aa_ref == aa_alt :
                        aa_change = "Syn"

                    domain = covid.get_s_domain(int(aa_position))

                    spike_mutations.append([nt_position,nt_ref+">"+nt_alt,aa_change,domain])

        df_spike = pd.DataFrame(spike_mutations)
        df_spike.columns = ["Genomic_locus","Nucleotide_change","Amino-acid_change","Domain"]

        runid_colname = runid+"("+str(df_run.shape[0])+")"
        df_spike_per_run = df_spike.groupby(list(df_spike.columns)).size().rename(runid_colname).reset_index()

        per_run_dfs.append(df_spike_per_run)

    df_final = per_run_dfs[0]
    for dfr in per_run_dfs[1:]:
        df_final = pd.merge(df_final,dfr,on=['Genomic_locus', 'Nucleotide_change', 'Amino-acid_change', 'Domain'],how='outer')

    df_final['Genomic_locus'] = df_final['Genomic_locus'].astype(int)

    df_final.to_excel(out_dir+"spike_protein_snp_table_"+tag+"_.xlsx",index=False)
