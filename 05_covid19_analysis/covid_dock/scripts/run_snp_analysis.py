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


def run_alignment_analysis(params):
    
    log_msg("processing version 1 data..")
    #### get data from start to run 10
    df_start_10 = run_version_1_analysis()

    log_msg("processing version 2 data..")
    #### get data for run 11 onwards
    nta_dir = params["container"]+params["nta_input"]
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


def run():
    params = read_run_params()

    tag=params["run_tag"]
    out_dir = params["container"]+params["nta_output"]


    if params["source"] == "alignment":

        df = run_alignment_analysis(params)
        variant_screener.genomewide_mutations(df,out_dir,tag)
    
    df_muttable = pd.read_excel(out_dir+"1_genomic_mutations_"+tag+".xlsx")
    df_muttable.nucleotide_change = df_muttable.nucleotide_change.apply(ast.literal_eval)
    df_muttable.nt_aa_pair = df_muttable.nt_aa_pair.apply(ast.literal_eval)

    ###generate summary    
    variant_screener.variant_screener_summary(df_muttable,out_dir,tag)


    db= params["container"]+params["database"]+params["run_date"]

    ##### generate per variant summary
    variant_screener.variant_screener_per_variant_summary(df_muttable,db,out_dir,tag)

    ## generate_snp_table
    log_msg("generating snp table for s protein..")
    mutation.generate_snp_table(df_muttable,db,out_dir,tag)


                    

if __name__ == '__main__':

    run()