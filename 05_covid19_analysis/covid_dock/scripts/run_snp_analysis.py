import pandas as pd
import os
import sys
import ast
sys.path.append(os.path.dirname(__file__))
from gen_utils.gen_io import read_run_params,log_msg
from genomes_in_df import mutation 
import variant_screener


def run():
    params = read_run_params()
    run = params["current_run"]
    out_home = params["container"]+"output/" 
    out_dir = out_home+run+"/"

    ##### preprare for variant calling 
    df_muttable = pd.read_excel(out_dir+"1_genomic_mutations_"+run+".xlsx")

    if params["source"] == "alignment_current":
        df_muttable_prev = pd.read_excel(out_home+"/previous_genomic_table/1_genomic_mutations_all.xlsx")
        df_muttable = df_muttable.append(df_muttable_prev,ignore_index=True)
        df_muttable.drop_duplicates("strain",inplace=True)
        df_muttable.to_excel(out_dir+"1_genomic_mutations_all.xlsx",index=False)
        df_muttable.to_excel(out_home+"/previous_genomic_table/1_genomic_mutations_all_"+run+".xlsx",index=False)

    df_muttable.nucleotide_change = df_muttable.nucleotide_change.apply(ast.literal_eval)
    df_muttable.nt_aa_pair = df_muttable.nt_aa_pair.apply(ast.literal_eval)

    # #### generate per variant summary with pangolin
    db_pangolin = params["container"]+"pangolin_analysis/pangolin_lineage_assignment_for_pipeline_"+run+".xlsx"
    db = out_dir+"/database/5_strains_run_group_analysis_table_with_samplesheet.csv"

    variant_screener.variant_screener_per_variant_summary_with_pangolin(df_muttable,db,db_pangolin,out_dir,run)

    print("check...")
    print(df_muttable.head())
    # generate_snp_table
    log_msg("generating snp table for s protein..")
    df_muttable.columns=["strain","total_mutations","nucleotide_change","nt_aa_pair","deletions"]
    mutation.generate_snp_table(df_muttable,db,out_dir,run)

    variant_screener.variant_screener_summary_pangolin(df_muttable,out_dir,run)

if __name__ == '__main__':

    run()