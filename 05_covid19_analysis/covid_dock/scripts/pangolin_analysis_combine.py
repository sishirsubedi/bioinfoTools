import pandas as pd
from gen_utils.gen_io import read_run_params,log_msg

params = read_run_params()
home_dir = params["container"]+"pangolin_analysis/"
runs=params["runs"]
current_run = params["current_run"]

df = pd.DataFrame()
for run in runs:
    print("processing..."+run)
    df_cur = pd.read_csv(home_dir+run+"/lineage_report_"+run+"_v3_1.csv")
    df_cur = df_cur[['taxon', 'lineage']]
    df_cur["consensus_run"] = run
    df = df.append(df_cur,ignore_index=True)

df["taxon"] = [x.replace("Consensus_","").replace(".ivar_threshold_0.6_quality_20_|_One_Codex_consensus_sequence","") for x in df["taxon"] ]

df["taxon"] = [x.replace("hCoV-19/USA/TX-HMH-","").replace("/2020","") for x in df["taxon"] ]

df.drop_duplicates("taxon",keep="last",inplace=True)

df = df[df["taxon"].str.contains("MCoV")]

df[['taxon','lineage','consensus_run']].to_excel(home_dir+"pangolin_lineage_assignment_for_pipeline_"+current_run+".xlsx",index=False)
