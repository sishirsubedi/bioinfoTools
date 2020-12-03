import pandas as pd
import alignment_analysis

########################################################

###read orders########
ref_dir="/home/tmhsxs240/COVID_19/reference/"
df_db = pd.read_csv(ref_dir+"1_orders.csv",dtype = {'MRN': str,'ORDER_ID': str},low_memory=False)
df_db.drop_duplicates("ORDER_ID",inplace=True)
df_db = df_db[['MRN','ORDER_ID','COLLECTION_DT']]


###read log files#####
df= pd.read_excel(ref_dir+"3_CoVID Sample Log.xlsx")
df = df[['Musser Lab No.', 'Full Order Number.2']]
df.columns = ['MCoVNumber', 'ORDER_ID']

dfjoin = pd.merge(df,df_db, on="ORDER_ID",how="left",indicator=True)
dfjoin._merge.value_counts()


####need to adress missing these in orders.csv 
####these are all training runs so ok to ignore these
#        MCoVNumber    ORDER_ID  MRN COLLECTION_DT     _merge
# 742      MCoV-743  I209005171  NaN           NaN  left_only
# 1218    MCoV-1219  I221003705  NaN           NaN  left_only
# 12122  MCoV-12123  I716012210  NaN           NaN  left_only
# 12290  MCoV-12291  I716012193  NaN           NaN  left_only
# 12787  MCoV-12788  I619009616  NaN           NaN  left_only
# 12823  MCoV-12824  I619009638  NaN           NaN  left_only
# 12866  MCoV-12867  I619009651  NaN           NaN  left_only

dfjoin = dfjoin[dfjoin._merge=="both"]
dfjoin = dfjoin[[ 'MCoVNumber','ORDER_ID','MRN','COLLECTION_DT']]
dfjoin.columns =['Strain', 'ORDER_ID','MRN','COLLECTION_DT']


dfjoin.COLLECTION_DT = pd.to_datetime(dfjoin.COLLECTION_DT)

def assign_group(cdate):
    cdate = pd.Timestamp(cdate)
    gr = ""
    if cdate < pd.Timestamp(2020,5,12):
        gr="wave_1"
    elif cdate > pd.Timestamp(2020,5,12) and cdate < pd.Timestamp(2020,8,15):
        gr="wave_2"
    elif cdate > pd.Timestamp(2020,8,15):
        gr="trough"
    return gr


dfjoin["AnalysisGroup"] = list(map(assign_group,dfjoin.COLLECTION_DT.values))
dfjoin.to_excel(ref_dir+"4_1_Curated_MCOV_MRN_Strains.xlsx",index=False)

#### read alignment ######
df_alignment = alignment_analysis.getStrains("/home/tmhsxs240/COVID_19/data/10_23/Houston.Oct.clean.fa")


#######################################

dfjoin_selected = pd.merge(df_alignment,dfjoin,how="left",on="Strain",indicator=True)

## remove
dfjoin_selected = dfjoin_selected[dfjoin_selected._merge=="both"]
dfjoin_selected = dfjoin_selected[dfjoin_selected.Strain != "MCoV-1255"]
dfjoin_selected = dfjoin_selected[dfjoin_selected.Strain != "MCoV-1343"]
# dfjoin_selected = dfjoin_selected[dfjoin_selected.Strain != "MCoV-1483"]

dfjoin_selected.to_csv(ref_dir+"4_2_Curated_MCOV_MRN_Strains_alignment_strains.csv",index=False)

################################################################################

#### generate strain count per patient 
t = dfjoin_selected.groupby(["MRN"])["Strain"].count().reset_index()
t2 = dfjoin_selected.groupby(["MRN"])["Strain"].agg({"Strain":'/'.join})
t2.reset_index(inplace=True)
t2["Strain Count"] = t["Strain"]
t2.to_excel(ref_dir+"4_3_Patients_strain_count_from_NTA.xlsx",index=False)



alignment_analysis.tableGenomicMutations("/home/tmhsxs240/COVID_19/data/10_23/Houston.Oct.clean.fa",ref_dir+"5_StrainsAndMutations_entire_genome.csv")
df_mutations_per_strain = pd.read_csv(ref_dir+"5_StrainsAndMutations_entire_genome.csv")


df_curated = pd.read_csv(ref_dir+"4_Curated_MCOV_MRN_Strains_5085.csv")
df_curated.columns =['Strain', 'ORDER_ID','MRN','COLLECTION_DT']
df_curated['COLLECTION_DT'] = pd.to_datetime(df_curated['COLLECTION_DT'])


dfjoin_selected = pd.merge(df_mutations_per_strain,df_curated,how="left",on="Strain",indicator=True)
dfjoin_selected = dfjoin_selected[dfjoin_selected._merge=="both"]

df_mrn_all = pd.DataFrame()
for mrn in dfjoin_selected.MRN.unique():
    df_mrn = dfjoin_selected[dfjoin_selected.MRN==mrn]
    df_mrn = df_mrn.sort_values('COLLECTION_DT')
    collection_dates = df_mrn.COLLECTION_DT.values
    collection_count = len(collection_dates)

    if collection_count >1:
        df_mrn_all = df_mrn_all.append(df_mrn)

df_mrn_all = df_mrn_all[['MRN','ORDER_ID','COLLECTION_DT','Strain', 'mutation_count', 'mutation_info' ]]


df_mrn_all["mutation_count"] = [str(x) for x in df_mrn_all.mutation_count]
t = df_mrn_all.groupby(["MRN"])["Strain"].count().reset_index()
t2 = df_mrn_all.groupby("MRN",as_index=False).agg({"Strain":'/'.join,"mutation_count":'/'.join})
t2.reset_index(inplace=True)
t2["Strain Count"] = t["Strain"]
t2 = t2[['MRN', 'Strain Count','mutation_count', 'Strain']]

#export main file
writer = pd.ExcelWriter("6_multiple_strains_snp_analysis.xlsx", engine='xlsxwriter')
df_mrn_all.to_excel(writer, sheet_name='StrainsInRow', index=False)
t2.to_excel(writer, sheet_name='MRNinRow', index=False)
writer.save()
