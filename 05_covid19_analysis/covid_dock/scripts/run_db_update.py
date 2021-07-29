import pandas as pd
import os
import sys
sys.path.append(os.path.dirname(__file__))
from gen_utils.gen_io import read_run_params,log_msg,bashCommunicator
from alignment.alignment_tools import strainFromAlignment
from gen_utils.gen_io import read_run_params
from genomes_in_df.filter import filterMainAlignment


def generate_mcov_db(ref_dir):
    
    orders =  "orders.csv"
    patients = "patients.csv"
    sample_log = "CoVID_Sample_Log.xlsx"
    out_file = "4_Curated_MCOV_MRN_Strains.xlsx"
    #######################
    bashCommunicator( "wget -P "+ ref_dir +" ")
    df_orders = pd.read_csv(ref_dir+orders, dtype = {'MRN': str,'ORDER_ID': str},low_memory=False)
    df_orders.drop_duplicates("ORDER_ID",inplace=True)
    df_orders = df_orders[['MRN','ORDER_ID','COLLECTION_DT','ORDERING_CLINIC_ID',
    'ORDERING_CLINIC_NAME', 'FACILITY','ADMISSION_DT','DISCHARGE_DT','HIS_PATIENT_TYPE']]

    bashCommunicator( "wget -P "+ ref_dir +" ")
    df_patients = pd.read_csv(ref_dir+patients, dtype = {'MRN': str,'ORDER.ORDER_ID': str},low_memory=False)
    df_patients = df_patients[['MRN', 'ZIP', 'ORDER.ORDER_ID']]
    df_patients.rename(columns={'ORDER.ORDER_ID':'ORDER_ID'},inplace=True)
    df_patients.drop_duplicates("ORDER_ID",inplace=True)
    df_patients['ZIP'] = [str(x)[0:5] for x in df_patients['ZIP']] 

    bashCommunicator("cp "+params["container"]+"temp/"+sample_log +" "+ ref_dir)
    df= pd.read_excel(ref_dir+sample_log)
    df = df[['Musser Lab No.', 'Full Order Number.2']]
    df.columns = ['MCoVNumber', 'ORDER_ID']

    dfjoin = pd.merge(df_orders,df_patients[['MRN','ZIP']],on='MRN',how='left')
    dfjoin = pd.merge(df,dfjoin,on='ORDER_ID',how='left',indicator=True)
    dfjoin.rename(columns={'_merge':'sequence_order_patient_data_match'},inplace=True)
    dfjoin.to_excel(ref_dir+out_file,index=False)


def mcov_order_id_dup_check(ref_dir):

    df = pd.read_excel(ref_dir+"4_Curated_MCOV_MRN_Strains.xlsx")

    df_dups = df[df["ORDER_ID"].isin(\
                df.groupby("ORDER_ID")\
                .agg('size')\
                .rename("mcov_count")\
                .reset_index()\
                .query("mcov_count>1")["ORDER_ID"])]\
                .sort_values("ORDER_ID")
    
    df_dups["dup_count"]=df_dups.groupby("ORDER_ID").transform('count')["MCoVNumber"]
    df_dups.to_excel(ref_dir+"4_qc_duplicate_order_ids.xlsx",index=False)


def generate_strain_ids(params,ref_dir):

    nta_input = params["container"]+"input_alignment/"

    all_strains = []
    for run in sorted(params["nta_group_v1"]):

        if "reference" in run:
            continue

        nta_file = params['nta_group_v1'][run].split(',')[0]
        all_strains.extend(strainFromAlignment(nta_input+nta_file,run))

    for run in sorted(params["nta_group_v2"]):

        if "reference" in run:
            continue
        
        nta_file = params['nta_group_v2'][run]
        all_strains.extend(strainFromAlignment(nta_input+nta_file,run))


    df = pd.DataFrame(all_strains)
    df = filterMainAlignment(df)
    df.columns = ["strain","run_group_analysis"]

    df_db = pd.read_excel(ref_dir+"4_Curated_MCOV_MRN_Strains.xlsx")

    df_db = df_db[df_db.sequence_order_patient_data_match=="both"]

    dfjoin = pd.merge(df_db,df,left_on="MCoVNumber",right_on="strain",how="outer",indicator=True)

    dfjoin.rename(columns={'_merge':'sequence_order_patient_alignment_data_match'},inplace=True)

    dfjoin = dfjoin[dfjoin.sequence_order_patient_alignment_data_match.isin(["both","left_only"])]

    dfjoin.loc[dfjoin.sequence_order_patient_alignment_data_match.isin(["left_only"]),["run_group_analysis"]]="pending_rungroup"

    dfjoin.to_csv(ref_dir+"5_strains_run_group_analysis_table.csv",index=False)


def update_run_id(x,y,z):
    
    new_id=""

    if x=="2_published":
        new_id="Run00_1_published"
    elif x=="3_after-published-run-1-4" and pd.isnull(y):
        new_id="Run00_2_afterpub_ONT"
    elif x=="pending_rungroup" and pd.notnull(y):
        new_id= y+"_low_quality"
    elif pd.isnull(y):
        new_id= "low_quality"
    elif x=="pending_rungroup" and pd.isnull(y):

        if z<pd.to_datetime("2020-12-30"):
            new_id = "low_quality"
        else:
            new_id = x    
    else:
        new_id=y
    
    return new_id


def update_all_LQ_in_run(dfjoin):
    count_table = dfjoin.run_id_seq.value_counts()
    check_runs = list(count_table[count_table>700].index)

    pending = []
    for run in check_runs:
        if "Run" in run and "low_quality" in run:
            pending.append(run)

    if len(pending) > 0:
        for run in pending:
            dfjoin.loc[dfjoin.run_id_seq==run,["run_id_seq"]]= run.replace("low_quality","")+"pending_analysis"
    
    return dfjoin


def add_quality(x):
    if "low_quality" in x:
        return "LQ"
    elif "pending" in x:
        return "Pending"
    else:
        return "HQ"


def generate_rungroup_from_samplesheet(ref_dir):

    from pathlib import Path
    path = Path(params["container"]+"input_samplesheet/")

    df_combine = pd.DataFrame()
    for samplesheet in path.rglob('*.csv'):
        run = samplesheet.name.split('_')[0]+"_"+samplesheet.name.split('_')[1]
        df =  pd.read_csv(samplesheet,skiprows=20,encoding= 'unicode_escape')
        df = df[["Sample_ID"]]
        if len(run.split('_')[1])<2:
            run = run.replace("run_","run_0")
        run = run.replace("r","R")
        df["run_id"]=run
        df_combine = df_combine.append(df)

    df_db = pd.read_csv(ref_dir+"5_strains_run_group_analysis_table.csv")

    dfjoin = pd.merge(df_db,df_combine,left_on="MCoVNumber",right_on="Sample_ID",how="outer",indicator=True)

    dfjoin.rename(columns={'_merge':'sequence_to_alignment_samplesheet_data_match'},inplace=True)

    dfjoin = dfjoin[dfjoin.sequence_to_alignment_samplesheet_data_match.isin(["both","left_only"])]
        
    dfjoin.COLLECTION_DT = pd.to_datetime(dfjoin.COLLECTION_DT)

    dfjoin["run_id_seq"] = list(map(update_run_id,dfjoin.run_group_analysis,dfjoin.run_id,dfjoin.COLLECTION_DT))

    dfjoin = update_all_LQ_in_run(dfjoin)

    dfjoin["quality"] = list(map(add_quality,dfjoin.run_id_seq))

    dfjoin["run_id_final"] = [x.replace("_low_quality","") for x in dfjoin["run_id_seq"] ]

    dfjoin.groupby(["run_group_analysis","run_id_seq"]).agg('count')["MCoVNumber"].reset_index().sort_values("run_id_seq").to_excel(ref_dir+"run_samplesheet_match_.xlsx",index=False)
    
    dfjoin.sort_values("run_id",inplace=True)

    dfjoin.drop_duplicates("MCoVNumber",keep="last",inplace=True)

    dfjoin.to_csv(ref_dir+"5_strains_run_group_analysis_table_with_samplesheet.csv",index=False)

        
params = read_run_params()
run = params["current_run"]
bashCommunicator("mkdir -p " +params["container"]+"output/"+run+"/database/")
bashCommunicator("cp "+params["container"]+"temp/"+run+".wuhan1.fa "+ params["container"]+"input_alignment/ " )
bashCommunicator("cp "+params["container"]+"temp/"+run+"_set1_SampleSheet.csv "+ params["container"]+"input_samplesheet/ " )
ref_dir = params["container"]+"output/"+run+"/database/"

generate_mcov_db(ref_dir)
generate_strain_ids(params,ref_dir)
generate_rungroup_from_samplesheet(ref_dir)