import pandas as pd
import alignment_analysis

########################################################

###read orders########
ref_dir="/home/tmhsxs240/COVID_19/reference/"
df_db = pd.read_csv(ref_dir+"1_orders.csv")
df_db.drop_duplicates("ORDER_ID",inplace=True)
df_db = df_db[['MRN','ORDER_ID','COLLECTION_DT']]
df_db.MRN = df_db.MRN.astype(str)
df_db.ORDER_ID = df_db.ORDER_ID.astype(str)


###read log files#####
df= pd.read_excel(ref_dir+"3_CoVID Sample Log.xlsx")
df = df[['Musser Lab No.', 'Full Order Number.2']]
df.columns = ['MCoVNumber', 'ORDER_ID']
df.MCoVNumber = df.MCoVNumber.astype(str)
df.ORDER_ID = df.ORDER_ID.astype(str)

dfjoin = pd.merge(df,df_db, on="ORDER_ID",how="left",indicator=True)
dfjoin._merge.value_counts()
dfjoin = dfjoin[dfjoin._merge=="both"]
dfjoin = dfjoin[[ 'MCoVNumber','ORDER_ID','MRN','COLLECTION_DT']]
dfjoin.columns =['Strain', 'ORDER_ID','MRN','COLLECTION_DT']
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
