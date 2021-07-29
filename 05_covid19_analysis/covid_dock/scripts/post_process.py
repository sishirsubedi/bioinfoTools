import pandas as pd
import ast
import sys
from gen_utils.gen_io import read_run_params,log_msg
 
def filter_spike_mutation_strainwise(df):
    strains=[]
    for indx,row in df.iterrows():
        spike =[]
        if pd.notnull(row.nt_aa_pair):
            for mut_pair in ast.literal_eval(row.nt_aa_pair):
                if ":S-" in mut_pair:
                    protein=mut_pair.split(':')[1]
                    aa_change=protein.split('-')[1]
                    if aa_change[0] != aa_change[len(aa_change)-1]:
                        spike.append(protein)
        strains.append(spike)
    return strains


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

    dfvar['spike-dels-from-var']=strains

    dfjoin = pd.merge(df,dfvar[['strain','spike-dels-from-var']],how="left",left_on="MCoVNumber",right_on="strain")
    return dfjoin


def curate_snvs_due_to_deletion_from_var(df):

    strains=[]
    for indx,row in df.iterrows():
        if isinstance(row['spike-dels-from-var'],float):
            strains.append("")
        elif "F157-" in row['spike-dels-from-var']:
            strains.append("E156G")
        else:
            strains.append("")
    df['spike-snvs-due-to-deletion']=strains

    return df


def analyze_voi(GENOMIC_FILE,META_FILE,VAR_FILE,out_dir,TAG):
    df = pd.read_csv(META_FILE)

    dfdb = pd.read_excel(GENOMIC_FILE)
    dfjoin = pd.merge(df,dfdb,how="left",left_on="MCoVNumber",right_on="strain")

    dfjoin["spike-snvs-from-alignment"] = filter_spike_mutation_strainwise(dfjoin)
    dfjoin= filter_deletion_from_var(dfjoin,VAR_FILE)
    dfjoin = curate_snvs_due_to_deletion_from_var(dfjoin)

    dfjoin = dfjoin[dfjoin.quality=="HQ"]

    dfjoin = dfjoin[['MCoVNumber', 'ORDER_ID', 'MRN', 'COLLECTION_DT', 'ZIP', 
       'run_id_final',  'variant', 'spike-snvs-from-alignment',
       'spike-dels-from-var', 'spike-snvs-due-to-deletion']]

    for col in ['spike-snvs-from-alignment','spike-dels-from-var', 'spike-snvs-due-to-deletion']:
        dfjoin[col] = [','.join(i) if isinstance(i, list) else i for i in dfjoin[col]]

    dfjoin.to_excel(out_dir+TAG+"_source_file.xlsx",index=False)


#############################################

params = read_run_params()
run = params["current_run"]
out_home = params["container"]+"output/" 
out_dir = out_home+run+"/"

GENOMIC_FILE = out_dir+"1_genomic_mutations_all.xlsx"
META_FILE = out_dir+"4_mcov_strain_variant_map_covid_pangolin_db_input_"+run+".csv"
VAR_FILE= params["container"]+"variant_analysis/final_result_with_var_assignment_"+run+".xlsx"


############## combine #######################################
analyze_voi(GENOMIC_FILE,META_FILE,VAR_FILE,out_dir,run)

