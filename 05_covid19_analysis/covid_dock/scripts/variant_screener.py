import pandas as pd
from alignment.alignment_tools import analyzeEntireGenomewithReference
from gen_utils.gen_io import read_run_params,log_msg
from genomes_in_df.filter import filterMainAlignment
from genomes_in_df import mutation 
import gen_utils.gen_covid as covid
import xlsxwriter as xl

def generate_df_from_multiple_alignments(alignment_files):

    df = pd.DataFrame()
    for alignment_file in alignment_files:
        df_af = analyzeEntireGenomewithReference(alignment_file)
        df_af = filterMainAlignment(df_af)
        df = df.append(df_af,ignore_index=True)
        print(df.shape)

    df = filterMainAlignment(df)
    
    log_msg("completed processing all the version 2 alignment files..")
    print(df.head())

    return df

def genomewide_mutations(df,out_dir,tag,deletion=True):

    log_msg("generating genome-wide mutations..")

    df_muttable = mutation.genomeWideMutationTable(df)

    df_muttable = mutation.add_aa_mutations_column(df_muttable,out_dir,tag)

    if deletion:
        log_msg("generating genome-wide Deletions..")

        df_muttable = mutation.add_aa_deletions_column(df_muttable,out_dir,tag)

    df_muttable.to_excel(out_dir+"1_genomic_mutations_"+tag+".xlsx",index=False)


def genomewide_mutations_outside(df,out_dir,tag,out_tag):

    log_msg("generating genome-wide mutations..for outside samples")

    df_out = df[df['id'].str[0:4].isin([out_tag,'MN90'])] 

    df_muttable = mutation.genomeWideMutationTable(df_out)

    df_muttable = mutation.add_aa_mutations_column(df_muttable,out_dir,tag)

    df_muttable.to_excel(out_dir+"1_genomic_mutations_ocov_"+tag+".xlsx",index=False)


def get_variants_match_mutations(df,variants):
    all_variants=[]
    for indx,row in df.iterrows():
        sel_variants = [x.split(':')[1] for x in row["nt_aa_pair"] if x.split(':')[1] in variants.values()]
        if len(sel_variants) >0:
            all_variants.append( sel_variants)
        else:
            all_variants.append([])
    return [len(x) for x in all_variants ],[",".join(x) for x in all_variants]

def variant_screener_summary(df_muttable,out_dir,tag):

    log_msg("screening variants..")

    var_dict = covid.get_covid_variants()
    
    for variant in var_dict:
        total_snvs = str(len(var_dict[variant]))    
        df_muttable[variant] ,df_muttable[variant+"(SNVs:"+total_snvs+")"] =\
        get_variants_match_mutations(df_muttable,var_dict[variant])
    
    df_muttable.to_excel(out_dir+"2_variants_screen_summary_"+tag+".xlsx",index=False)

def variant_screener_summary_pangolin(df_muttable,out_dir,tag):

    log_msg("screening variants..")

    var_dict = covid.get_covid_variants()
    
    for variant in var_dict:
        total_snvs = str(len(var_dict[variant]))    
        df_muttable[variant] ,df_muttable[variant+"(SNVs:"+total_snvs+")"] =\
        get_variants_match_mutations(df_muttable,var_dict[variant])
    
    df_muttable.to_excel(out_dir+"2_variants_screen_summary_"+tag+".xlsx",index=False)

def write_excel_table(writer,sheet_name):

    workbook  = writer.book
    my_format = workbook.add_format()


    worksheet = writer.sheets[sheet_name]

    my_format.set_align('center')

    if sheet_name == "Summary":
        worksheet.set_column('A:M', 22,my_format)
    else:
        worksheet.set_column('A:H', 20,my_format)

def variant_screener_per_variant_summary_with_pangolin(df,db,db_pangolin,out_dir,tag):

    log_msg("generating per variant information..")

    df_db = pd.read_csv(db)
    df.rename(columns={"strain":"MCoVNumber"},inplace=True)
    dfjoin = pd.merge(df,df_db,on="MCoVNumber",how="left",indicator=True)
    dfjoin.rename(columns={'_merge':'variant_db_Match'},inplace=True)

    dfjoin = dfjoin[dfjoin.variant_db_Match=="both"]
    dfjoin.COLLECTION_DT = pd.to_datetime(dfjoin.COLLECTION_DT)

    ### add pangolin
    df_pangolin = pd.read_excel(db_pangolin)
    dfjoin = pd.merge(dfjoin,df_pangolin,right_on="taxon",left_on="MCoVNumber",how="left")
    print(dfjoin.head())


    writer = pd.ExcelWriter(out_dir+"3_variant_screening_per_variant_summary_pangolin_"+tag+".xlsx", engine='xlsxwriter',
    datetime_format='yyyy-mm-dd hh:mm:ss', date_format='yyyy-mm-dd')

    vois = covid.covid_variants.pangolin

    summary_variants = {}
    mcov_variant = []
    for variant in vois:
    
        df_variant = dfjoin[dfjoin["lineage"]==variant]

        df_variant = df_variant[['MCoVNumber', 'ORDER_ID', 'MRN','COLLECTION_DT', 
            'ORDERING_CLINIC_ID','ORDERING_CLINIC_NAME',
            'FACILITY','ADMISSION_DT','DISCHARGE_DT','HIS_PATIENT_TYPE',
            'ZIP','lineage','run_id_seq']]

        df_variant.sort_values("run_id_seq",inplace=True)

        for indx,row in df_variant.iterrows():
            mcov_variant.append(
                [
                    row['MCoVNumber'], variant
                ]
            )

        df_variant.to_excel(writer, sheet_name=variant,index=False)

        write_excel_table(writer,variant)

        summary_variants[variant] = df_variant.groupby("run_id_seq").size().to_dict()
    
    ###save mcov variants file for covid database
    df_mcov_variant = df_pangolin[["taxon","lineage"]]
    df_mcov_variant.columns = ['strain','variant']
    dfjoin_mcov_variant = pd.merge(df_db,df_mcov_variant,left_on="MCoVNumber",right_on="strain",how="left")
    
    
    mcov_first_variant_per_patient = dfjoin_mcov_variant[dfjoin_mcov_variant.quality=="HQ"].sort_values(by='COLLECTION_DT', ascending=True).drop_duplicates(['variant', 'MRN'])['MCoVNumber'].values
    dfjoin_mcov_variant['Is_First_For_Patient'] = [1 if x in mcov_first_variant_per_patient else 0 for x in dfjoin_mcov_variant["MCoVNumber"]]    
    dfjoin_mcov_variant.to_csv(out_dir+"4_mcov_strain_variant_map_covid_pangolin_db_input_"+tag+".csv",index=False)

    sel_runs = sorted([x for x in dfjoin.run_id_seq.unique() if "low_quality" not in x])
    
    summary_table = []
    summary_column = ["Run","Total","Start-Date","End-Date"]
    for variant in sorted(summary_variants):
        summary_column.append(variant)
    summary_table.append(summary_column)

    for runid in sel_runs:
        
        run_row = []

        df_runid = dfjoin[dfjoin.run_id_seq==runid]
        total = df_runid.shape[0]
        start = df_runid.COLLECTION_DT.min().date()
        end = df_runid.COLLECTION_DT.max().date()

        run_row.append(runid)
        run_row.append(total)
        run_row.append(start)
        run_row.append(end)

        for variant in sorted(summary_variants):
            try:
                run_row.append(summary_variants[variant][runid])
            except:
                run_row.append("")

        summary_table.append(run_row)
    
    df_summary = pd.DataFrame(summary_table)
    df_summary.columns = df_summary.iloc[0,:]
    df_summary.drop(0, inplace=True)
    df_summary.to_excel(writer, sheet_name="Summary",index=False)
    write_excel_table(writer,"Summary")
    
    writer.close()
