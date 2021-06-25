import pandas as pd
import ast

def filter_mutation_of_interest(df,moi):
    keep=[]
    for indx,row in df.iterrows():
        for mut_pair in ast.literal_eval(row.nt_aa_pair):
            if moi in mut_pair:
                keep.append([row["strain"],mut_pair])
                break
    return keep
 
def filter_spike_mutation_strainwise(df):
    strains=[]
    for indx,row in df.iterrows():
        spike =[]
        for mut_pair in ast.literal_eval(row.nt_aa_pair):
            if ":S-" in mut_pair:
                protein=mut_pair.split(':')[1]
                aa_change=protein.split('-')[1]
                if aa_change[0] != aa_change[len(aa_change)-1]:
                    spike.append(protein)
        strains.append(spike)
    return strains

def filter_spike_mutation_haplotype(df):
    smut_dict = {}
    for indx,row in df.iterrows():
        pattern=""
        for mut_ in row["spike-snps"]:
            pattern +=mut_.split("-")[1]+"-"
        if pattern not in smut_dict:
            smut_dict[pattern] = []
            smut_dict[pattern].append(str(1))
            smut_dict[pattern].append(row["MCoVNumber"]+",")
            smut_dict[pattern].append(row["run_id_final"]+",")
        else:
            upd_val = int(smut_dict[pattern][0])+1
            smut_dict[pattern][0] = str(upd_val)
            smut_dict[pattern][1] += row["MCoVNumber"]+","
            smut_dict[pattern][2] += row["run_id_final"]+","
    return smut_dict


def filter_deletion_from_alignment(df):
    strains=[]
    for indx,row in df.iterrows():
        spike =[]
        for mut_pair in ast.literal_eval(row.nucleotide_change):
            if len(mut_pair)>=5 and "-" in mut_pair:
                if int(mut_pair[1:len(mut_pair)-1]) < 29675 :
                    spike.append(mut_pair)
        strains.append(spike)
    return strains

def filter_deletion_from_var(df,var_file):

    dfvar = pd.read_excel(var_file)

    strains=[]
    for indx,row in dfvar.iterrows():
        spike_del =[]
        for var in row['spike-var'].split(','):
            if "-" in var:
                spike_del.append(var)
        strains.append(spike_del)

    dfvar['spike-var-del']=strains

    dfjoin = pd.merge(df,dfvar[['strain','spike-var', 'spike-var-del']],how="left",left_on="MCoVNumber",right_on="strain")
    return dfjoin


def analyze_moi(GENOMIC_FILE,META_FILE,VAR_FILE,MOI):
    df = pd.read_excel(GENOMIC_FILE)
    keep = filter_mutation_of_interest(df,MOI)
    keep = [x[0] for x in keep]
    dfdb = pd.read_csv(META_FILE)
    dfdb = dfdb[dfdb.MCoVNumber.isin(keep)]

    ###
    dfdb = dfdb[dfdb.quality=="HQ"]
    dfdb.sort_values("MRN",inplace=True)
    dfdb.drop_duplicates("MRN",inplace=True)
    
    # #################optional - filter date
    # start_date = "1-1-2021"
    # end_date = "6-1-2021"
    # dfdb.COLLECTION_DT = pd.to_datetime(dfdb.COLLECTION_DT)
    # dfdb = dfdb[ dfdb.COLLECTION_DT>=pd.to_datetime(start_date)]
    # dfdb = dfdb[ dfdb.COLLECTION_DT<pd.to_datetime(end_date)]

    dfjoin = pd.merge(dfdb,df,how="left",left_on="MCoVNumber",right_on="strain")
    dfjoin["spike-snps"]=dfjoin["spike-snps"] = filter_spike_mutation(dfjoin)
    dfjoin["Interest"] = filter_mutation_of_interest_witn_mut(dfjoin,moi)
    dfjoin.to_excel(MOI+"_analysis.xlsx")


def analyze_voi(GENOMIC_FILE,META_FILE,VAR_FILE,VOI):
    df = pd.read_csv(META_FILE)

    df = df[df.quality=="HQ"]
    df = df[df.variant.notnull()]
    df = df[df.variant.str.contains(VOI)]

    dfdb = pd.read_excel(GENOMIC_FILE)
    dfjoin = pd.merge(df,dfdb,how="left",left_on="MCoVNumber",right_on="strain")

    dfjoin["spike-snps"] = filter_spike_mutation_strainwise(dfjoin)

    dfjoin= filter_deletion_from_var(dfjoin,VAR_FILE)

    # dfjoin["Interest"] = filter_mutation_of_interest_witn_mut(dfjoin,"S-E484")
    # dfjoin["deletions"] = filter_DEL_mutation(dfjoin)
    dfjoin.to_excel(VOI.replace(".","_")+"_analysis.xlsx")

    
    ##########get haplotpye
    res = filter_spike_mutation_haplotype(dfjoin)
    df_res = pd.DataFrame(res)
    df_res.to_excel(VOI.replace(".","_")+"_subvariant_analysis.xlsx")


#############################################
GENOMIC_FILE=""
META_FILE=""
VAR_FILE=""
MOI="S-E484"
VOI="B.1.617.2"
#############################################

############## analyze MOI#######################################
analyze_moi(GENOMIC_FILE,META_FILE,VAR_FILE,MOI)

############## analyze VOI#######################################
analyze_voi(GENOMIC_FILE,META_FILE,VAR_FILE,VOI)

