import pandas as pd
import os
import sys
sys.path.append(os.path.dirname(__file__))
from gen_utils.gen_io import read_run_params,log_msg
import process_version_1_data
import variant_screener


def run_version_1_analysis():
    return process_version_1_data.get_combine_df()


def run_alignment_analysis(params,mode):
    
    nta_dir = params["container"]+"input_alignment/"

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
        
        df_start_10 = df_start_10.append(df_11_18,ignore_index=True)

        df_start_10.rename(columns={0:"id"},inplace=True)
        print("finalized data")
        print(df_start_10.head())

        return df_start_10
    
    elif mode=="alignment_current":

        alignment_file = [nta_dir+params["alignment_current"]]

        df = variant_screener.generate_df_from_multiple_alignments(alignment_file)

        sel_indx = [0]
        for x in range(266,29675,1): sel_indx.append(x)
        df = df[sel_indx]
        df.rename(columns={0:"id"},inplace=True)
        print("finalized version 2 data")
        print(df.head())

        return df

def run():
    params = read_run_params()
    run = params["current_run"]
    out_home = params["container"]+"output/" 
    out_dir = out_home+run+"/"

    ##### generate genome wide variants
    if params["source"] == "alignment_all":

        df = run_alignment_analysis(params,"alignment_all")
        variant_screener.genomewide_mutations(df,out_dir,run)

    elif params["source"] == "alignment_current":

        df = run_alignment_analysis(params,"alignment_current")
        variant_screener.genomewide_mutations(df,out_dir,run)
                    

if __name__ == '__main__':

    run()