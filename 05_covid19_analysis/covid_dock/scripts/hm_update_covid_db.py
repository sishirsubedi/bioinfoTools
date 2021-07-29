import pandas as pd
import os
import sys
import ast
from gen_utils.gen_io import read_run_params,log_msg
import sqlalchemy
from sqlalchemy import create_engine, text
import pymysql

def getcon():
    params = read_run_params()
    SQLALCHEMY_DATABASE_URI = "mysql+pymysql://"+params["dbcred"]
    sqlEngine = sqlalchemy.create_engine(SQLALCHEMY_DATABASE_URI, echo=False)
    dbConnection = sqlEngine.connect()
    return dbConnection

def clearTable(tablename):
    dbConnection = getcon()
    dbConnection.execute(text("DELETE FROM " + tablename + " WHERE 1=1 "))
    dbConnection.close()

def saveToDatabase(df,tablename):
    dbConnection = getcon()
    df.to_sql(tablename, dbConnection, if_exists='append', index=False)
    dbConnection.close()

def upload_sample(sample_file):
    df = pd.read_csv(sample_file,dtype={
    'MCoVNumber': 'str',
    'ORDER_ID': 'str',
    'MRN': 'str'   
    }, parse_dates=['COLLECTION_DT'],low_memory=False)
    df = df[df['sequence_to_alignment_samplesheet_data_match']!='right_only']
    df = df[['MCoVNumber', 'ORDER_ID', 'MRN', 'COLLECTION_DT','ZIP', 'run_id_final','variant','quality','Is_First_For_Patient']]
    df.columns = ['mcov_id', 'order_id','mrn','collection_date','zipcode','run_group','variant','quality','is_first_for_patient']
    df.drop_duplicates('mcov_id',inplace=True)

    df = df[~df['mcov_id'].isnull()]
    # df = df[~df['rungroup'].isnull()]    
    saveToDatabase(df,"sample")

    return df

def upload_variant_growth(variant_growth_file):
    df = pd.read_csv(variant_growth_file, parse_dates=['collection_dt'],low_memory=False)
    saveToDatabase(df[['collection_dt','total_samples','variant_samples','variant']],"variant_growth")

def upload_mutations(df_sample,mutations_file):
    df_m = pd.read_excel(mutations_file)
    df_m = df_m[df_m.strain.isin(df_sample.mcov_id)]
    mutations = [] 
    for indx,row in df_m.iterrows():
        for mut in ast.literal_eval(row["nt_aa_pair"]):
            gene_aa = mut.split(':')[1]
            mutations.append([row["strain"],gene_aa.split('-')[0],gene_aa.split('-')[1]])
    
    df_mutations = pd.DataFrame(mutations,columns=['mcov_id','gene','mutation'])
    df_mutations.to_excel("sample_mutations.xlsx",index=False)
    saveToDatabase(df_mutations,"sample_mutations")

def load_variants(variant_name,variants):
    varlist =[]
    for variant in variants: 
        varlist.append([variant_name,variant.split('-')[0],variant.split('-')[1]])
    return varlist

def upload_variants():
    import gen_covid

    variants = gen_covid.get_covid_variants()

    varlist =[]
    for variant in variants:
        var_dict = variants[variant] 
        for snv in var_dict.values():
            varlist.append([variant,snv.split('-')[0],snv.split('-')[1]])
    
    df_variants = pd.DataFrame(varlist,columns=['variant_id','gene','mutation'])
    df_variants.to_excel("variant.xlsx",index=False)
    saveToDatabase(df_variants,"variant")

params = read_run_params()
run = params["current_run"]
out_dir = params["container"]+"output/"+run+"/"
dbcred = params["dbcred"]

sample_file = out_dir+"/4_mcov_strain_variant_map_covid_pangolin_db_input_"+run+".csv"
mutations_file = out_dir+"/1_genomic_mutations_all.xlsx"
variant_growth_file = out_dir+"/variant_growth.csv"

clearTable("sample_mutations")
clearTable("sample")
clearTable("variant_growth")

df_sample = upload_sample(sample_file)
upload_mutations(df_sample,mutations_file)
upload_variant_growth(variant_growth_file)
