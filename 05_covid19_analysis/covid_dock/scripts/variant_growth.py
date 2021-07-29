import pandas as pd
import numpy as np
from datetime import timedelta
import gen_utils.gen_covid as covid
from gen_utils.gen_io import read_run_params,log_msg
import sys

params = read_run_params()
run = params["current_run"]
out_dir = params["container"]+"output/"+run+"/"
df = pd.read_csv(out_dir+"4_mcov_strain_variant_map_covid_pangolin_db_input_"+run+".csv")
df = df[df.quality=="HQ"]
start_date = sys.argv[1]
end_date = sys.argv[2]

dfs=[]
for voi in covid.covid_variants.pangolin:
     print(voi)
     ###take unique patients with variant
     keep_mrns_variant = np.unique(df[df.variant==voi]["MRN"])
     df_mrns = df[df.MRN.isin(keep_mrns_variant)]
     df_mrns = df_mrns[df_mrns.variant==voi] ###important step--remove non b117 variant 
     df_mrns.sort_values("COLLECTION_DT",inplace=True)
     df_mrns.drop_duplicates("MRN",keep="first",inplace=True)

     keep_mrns_not_variant = np.unique(df[df.variant!=voi]["MRN"])
     df_mrns_not_variant = df[df.MRN.isin(keep_mrns_not_variant)]
     df_mrns_not_variant = df_mrns_not_variant[df_mrns_not_variant.variant!=voi]
     df_mrns_not_variant.sort_values("COLLECTION_DT",inplace=True)
     df_mrns_not_variant.drop_duplicates("MRN",keep="first",inplace=True)

     df_2 = df_mrns.append(df_mrns_not_variant,ignore_index=True)
     df_2.drop_duplicates("MRN",keep="first",inplace=True)

     df_2=df_2[['MCoVNumber','COLLECTION_DT','variant']]

     ########### add run 30 96 samples manually
     dfr30=pd.read_excel("/storage/scratch/covid/container/pangolin_analysis/run_30/96_samples_no_id/lineage_report_run_30_with_colletiondate.xlsx")
     df_2 = df_2.append(dfr30)

     df_2.drop_duplicates("MCoVNumber",keep="first",inplace=True)
     #####################################

     df_2.COLLECTION_DT = pd.to_datetime(df_2.COLLECTION_DT)
     df_2 = df_2[ (  (df_2.COLLECTION_DT>=pd.to_datetime(start_date)) &
               (df_2.COLLECTION_DT<=pd.to_datetime(end_date)) )
                ]
     df_2.sort_values("COLLECTION_DT",inplace=True)

     df_2.variant.fillna(0,inplace=True)
     #########################

     df_2.COLLECTION_DT = df_2.COLLECTION_DT.dt.date
     df_2.sort_values("COLLECTION_DT",inplace=True)
     df_2.variant = [1 if x==voi else 0 for x in df_2.variant]

     df_variant = df_2.groupby("COLLECTION_DT")["variant"].agg("sum").reset_index()
     df_count = df_2.groupby("COLLECTION_DT")["variant"].agg("count").reset_index()

     dates = pd.date_range(df_2.COLLECTION_DT.min(), (df_2.COLLECTION_DT.max() + timedelta(days=1) )-timedelta(days=1),freq='d')
     
     df_data = pd.DataFrame(dates)
     df_data.columns=["collection_dt"]
     df_data["date_step"]= [x for x in range(1,df_data.shape[0]+1,1)]
     df_data["total_samples"] = df_count.variant
     df_data["variant_samples"] = df_variant.variant
     df_data["variant_percent"]=[ (x/y)*100 for x,y in zip(df_data.variant_samples,df_data.total_samples)]
     df_data["variant"]=voi

     dfs.append(df_data)

df_final_data=dfs[0]
for df_data in dfs[1:]:
     df_final_data=df_final_data.append(df_data)

df_final_data.to_csv(out_dir+"/variant_growth.csv",index=False)

