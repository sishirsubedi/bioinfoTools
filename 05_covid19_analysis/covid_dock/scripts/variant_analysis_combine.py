import pandas as pd
from gen_utils.gen_io import read_run_params,log_msg,bashCommunicator

params = read_run_params()
home_dir = params["container"]+"variant_analysis/"
runs=params["runs"]
current_run=params["current_run"]

bashCommunicator("cp "+params["container"]+"temp/"+current_run+".var "+ home_dir+"runs/ " )

df = pd.DataFrame()
for run in runs:
    if run in ["run_old","run0_10","run11_15"]: 
        continue
    else:
        print("processing..."+run)
        df_cur = pd.read_csv(home_dir+"runs/"+run+".var",header=None,sep='\t')
        df_cur.columns = ['strain', 'num1','num2','spike-var']
        df_cur = df_cur[['strain','spike-var']]
        df_cur["run"] = run
        df = df.append(df_cur,ignore_index=True)

df.to_excel(home_dir+"final_result_with_var_assignment_"+current_run+".xlsx",index=False)
