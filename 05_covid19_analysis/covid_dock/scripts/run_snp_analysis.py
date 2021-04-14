import pandas as pd
import os
import sys
import ast
sys.path.append(os.path.dirname(__file__))
from gen_utils.gen_io import read_run_params,log_msg
import process_version_1_data
from genomes_in_df import mutation 
import variant_screener


def run_version_1_analysis():
    return process_version_1_data.get_combine_df()


def run_alignment_analysis(params,mode):
    
    nta_dir = params["container"]+params["nta_input"]

    if mode=="alignment_all":
        log_msg("processing version 1 data..")
        #### get data from start to run 10
        df_start_10 = run_version_1_analysis()

        log_msg("processing version 2 data..")
        #### get data for run 11 onwards
        nta_files = params["nta_group_v2"]

        alignment_files =[]
        log_msg("processing newer version alignments (after run11 novaseq - one codex pipeline) --")
        for nta_group in sorted(nta_files) :
            print(nta_group,nta_files[nta_group])
            alignment_files.append(nta_dir+nta_files[nta_group])

        df_11_18 = variant_screener.generate_df_from_multiple_alignments(alignment_files)

        sel_indx = [0]
        for x in range(266,29675,1): sel_indx.append(x)

        df_11_18 = df_11_18.loc[1:,sel_indx] ###remove reference
        df_11_18.rename(columns={0:"id"},inplace=True)
        print("finalized version 2 data")
        print(df_11_18.head())

        return df_start_10.append(df_11_18,ignore_index=True)
    
    elif mode=="alignment_current":

        alignment_file = [nta_dir+params["alignment_current"]]

        df = variant_screener.generate_df_from_multiple_alignments(alignment_file)

        sel_indx = [0]
        for x in range(266,29675,1): sel_indx.append(x)
        df.rename(columns={0:"id"},inplace=True)
        print("finalized version 2 data")
        print(df.head())

        return df

def run():
    params = read_run_params()

    tag=params["run_tag"]
    out_dir = params["container"]+params["nta_output"]

    ##### generate genome wide variants
    if params["source"] == "alignment_all":

        df = run_alignment_analysis(params,"alignment_all")
        variant_screener.genomewide_mutations(df,out_dir,tag)

    elif params["source"] == "alignment_current":

        df = run_alignment_analysis(params,"alignment_current")
        variant_screener.genomewide_mutations(df,out_dir,tag)


    ##### preprare for variant calling 
    df_muttable = pd.read_excel(out_dir+"1_genomic_mutations_"+tag+".xlsx")

    if params["source"] == "alignment_current":
        df_muttable_prev = pd.read_excel(out_dir+"/previous_genomic_table/1_genomic_mutations_all.xlsx")
        df_muttable = df_muttable.append(df_muttable_prev,ignore_index=True)
        df_muttable.drop_duplicates("strain",inplace=True)
        df_muttable.to_excel(out_dir+"1_genomic_mutations_all.xlsx",index=False)

    df_muttable.nucleotide_change = df_muttable.nucleotide_change.apply(ast.literal_eval)
    df_muttable.nt_aa_pair = df_muttable.nt_aa_pair.apply(ast.literal_eval)


    variant_screener.variant_screener_summary_pangolin(df_muttable,out_dir,tag)

    db= params["container"]+params["database"]+params["run_date"]

    # #### generate per variant summary with pangolin
    db_pangolin= params["container"]+params["pangolin"]
    variant_screener.variant_screener_per_variant_summary_with_pangolin(df_muttable,db,db_pangolin,out_dir,tag)

    # generate_snp_table
    log_msg("generating snp table for s protein..")
    df_muttable.columns=["strain","total_mutations","nucleotide_change","nt_aa_pair"]
    mutation.generate_snp_table(df_muttable,db,out_dir,tag)


                    

if __name__ == '__main__':

    run()