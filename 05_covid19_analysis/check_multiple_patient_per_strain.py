import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import AlignIO, SeqIO, Entrez
import ast
import pylab
import pandas as pd
import seaborn as sns
import numpy as np
import os
import re


df= pd.read_excel("3_MCoV Sample Log 6-22-20 for Paul.xlsx")
df = df[['Musser Lab No.', 'Full Order Number']]
df.columns = ['MCoVNumber', 'ORDER_ID']
df.MCoVNumber = df.MCoVNumber.astype(str)
df.ORDER_ID = df.ORDER_ID.astype(str)

df_db = pd.read_excel("1_covid_patient_orders.xlsx")
df_db2 = pd.read_excel("2_Surveillance.xlsx")
df_db = df_db.append(df_db2)
### keep only unique order ids
df_db.drop_duplicates("ORDER_ID",inplace=True)
df_db = df_db[['MRN','ORDER_ID','COLLECTION_DT']]
df_db.MRN = df_db.MRN.astype(str)
df_db.ORDER_ID = df_db.ORDER_ID.astype(str)

dfjoin = pd.merge(df,df_db, on="ORDER_ID",how="left",indicator=True)
dfjoin._merge.value_counts()
dfjoin = dfjoin[dfjoin._merge=="both"]
dfjoin = dfjoin[[ 'MCoVNumber','ORDER_ID','MRN','COLLECTION_DT']]
dfjoin.to_csv("4_Curated_MCOV_MRN_Strains.csv",index=False)



df_alignment = alignment_analysis.getStrains("/home/tmhsxs240/COVID_19/data/6_24/Houston.July1.clean--RedundantMRN.fa")

df_curated = pd.read_csv("4_Curated_MCOV_MRN_Strains.csv")
df_curated = df_curated[['MCoVNumber', 'ORDER_ID', 'MRN']]
df_curated.columns =['Strain', 'ORDER_ID','MRN']
dfjoin_selected = pd.merge(df_alignment,df_curated,how="left",on="Strain",indicator=True)
dfjoin_selected = dfjoin_selected[dfjoin_selected._merge=="both"]

t = dfjoin_selected.groupby(["MRN"])["Strain"].count().reset_index()
t2 = dfjoin_selected.groupby(["MRN"])["Strain"].agg({"Strain":'/'.join})
t2.reset_index(inplace=True)
t2["Strain Count"] = t["Strain"]
t2.to_excel("Patients_strain_count_from_NTA.xlsx",index=False)



alignment_analysis.tableGenomicMutations("/home/tmhsxs240/COVID_19/data/6_24/Houston.July1.clean--RedundantMRN.fa","5_StrainsAndMutations_entire_genome.csv")
df_mutations_per_strain = pd.read_csv("5_StrainsAndMutations_entire_genome.csv")


df_curated = pd.read_csv("4_Curated_MCOV_MRN_Strains.csv")
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

#####
##### number of changes and strains count
#
# from functools import partial
# def checkDiff(mlist,th):
#     mutation_count = mlist.split('/')
#     first_num = int(mutation_count[0])
#     diff = 0
#     mut_diff=[]
#     for n in mutation_count[1:]:
#         mut_diff.append(abs(first_num-int(n)))
#     if th==np.max(mut_diff):
#         diff = 1
#     return diff
#
# for i in range(0,20,1):
#     col_name = "mutation_diff_"+str(i)
#     t2[col_name] = list(map(partial(checkDiff,th=i),t2.mutation_count.values))
#
#
# sum_table = t2[["mutation_diff_0","mutation_diff_1","mutation_diff_2","mutation_diff_3","mutation_diff_4",
# "mutation_diff_5","mutation_diff_6","mutation_diff_7",
# "mutation_diff_8","mutation_diff_9","mutation_diff_10",
# "mutation_diff_11","mutation_diff_12","mutation_diff_13",
# "mutation_diff_14","mutation_diff_15"
# ]].sum()
#
# df_sum_table = pd.DataFrame(sum_table)
# df_sum_table.columns = ["Total_strains"]
# df_sum_table["percent"] =  round((df_sum_table["Total_strains"]/t2.shape[0])*100,2)


#### now look at patients with two strains with one mutation

# df_2strain_1mutation = t2[t2["Strain Count"]==2]

df_1mutation = t2.copy()

df_1mutation = df_1mutation[['MRN', 'Strain Count', 'mutation_count', 'Strain']]

df_1mutation_join = pd.merge(df_1mutation,df_mrn_all, on="MRN",how="left")

df_1mutation_join = df_1mutation_join[['MRN', 'Strain Count', 'mutation_count_x',
'Strain_x', 'ORDER_ID', 'Strain_y', 'mutation_count_y', 'mutation_info']]


mrn_changed_mutation ={}

for mrn in df_1mutation_join.MRN.unique():

    df_patient = df_1mutation_join[df_1mutation_join.MRN==mrn]

    all_mutations = []
    for indx,mut_list in df_patient.mutation_info.items():
        mut_list = ast.literal_eval(mut_list)
        for mut in mut_list:
            if mut not in all_mutations:
                all_mutations.append(mut)

    changed_mutations = []
    for indx,mut_list2 in df_patient.mutation_info.items():
        mut_list2 = ast.literal_eval(mut_list2)
        for mut2 in all_mutations:
            if mut2 not in mut_list2 and mut2 not in changed_mutations:
                changed_mutations.append(mut2)

    mrn_changed_mutation[mrn] = changed_mutations

df_mrn_changed_mutation = pd.DataFrame(list(mrn_changed_mutation.items()),columns = ['MRN','changed_mutations'])
df_mrn_changed_mutation["total_change"] = [len(x) for x in df_mrn_changed_mutation.changed_mutations]
df_mrn_changed_mutation = df_mrn_changed_mutation[['MRN','total_change','changed_mutations']]


df_mrn_changed_mutation_join = pd.merge(df_1mutation_join,df_mrn_changed_mutation, on="MRN",how="left")
df_mrn_changed_mutation_join.to_excel("7_multiple_strains_snp_analysis_changed_mutations_only.xlsx",index=False)


##########
strains=2
df2 = pd.DataFrame(df_mrn_changed_mutation_join[df_mrn_changed_mutation_join["Strain Count"]==strains].total_change.value_counts())
df2["patients"] = [int(x/strains) for x in df2.total_change]
df2.reset_index(inplace=True)
df2.sort_values("index",inplace=True
df2[["index","patients"]]


##### annotation variants

vcf_file = []
for indx,row in df_mrn_changed_mutation_join.iterrows():
    for mutation in row["changed_mutations"]:
        if mutation != "":
            vcf_file.append(["NC_045512.2",mutation[1:len(mutation)-1],".",mutation[0],mutation[len(mutation)-1],"100.0","PASS","INFO"])
df_vcf_file = pd.DataFrame(vcf_file)
df_vcf_file[1] = df_vcf_file[1].astype(int)
df_vcf_file = df_vcf_file.drop_duplicates()
df_vcf_file = df_vcf_file.sort_values(1)
df_vcf_file.to_csv("covid.vcf",sep='\t',header=None,index=False)


sudo java -Xmx4g -jar /opt/snpeff/snpeff_covid/snpEff/snpEff.jar -v   NC_045512.2   covid.vcf > covid.vcf.snpeff

#################################################
# count mutations based on number of strains
#################################################


VEP_FILE="covid.vcf.snpeff"
header=[]
with open(VEP_FILE) as myfile:
    header = [next(myfile) for x in range(3)]
columns_line=str(header[2].strip().split(':')[1]).split('|')

df_snpeff = pd.read_csv(VEP_FILE,sep='\t',skiprows=5,header=None)
df_snpeff.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
df_snpeff["HGMD_INFO"]= [x.split(";")[0] for x in df_snpeff.INFO]
df_snpeff["VEP_INFO"]= [x.split(";")[1].split(',')[0] for x in df_snpeff.INFO]
df_snpeff[columns_line] = df_snpeff['VEP_INFO'].str.split('|',expand=True)



writer = pd.ExcelWriter("Final_snp_analysis_strains_per_patient.xlsx", engine='xlsxwriter')

for strain_count in [2,3,4,5,6]:

    df_mrn_changed_mutation_join_all_tc_strain_count_strain = df_mrn_changed_mutation_join[df_mrn_changed_mutation_join["Strain Count"]==strain_count]

    df_mrn_changed_mutation_join_all_tc_strain_count_strain = df_mrn_changed_mutation_join_all_tc_strain_count_strain[df_mrn_changed_mutation_join_all_tc_strain_count_strain["total_change"]>0]

    df_mrn_changed_mutation_join_all_tc_strain_count_strain = df_mrn_changed_mutation_join_all_tc_strain_count_strain.drop_duplicates(['MRN'])


    new_mutations = []

    for indx,row in df_mrn_changed_mutation_join_all_tc_strain_count_strain.iterrows():
        mrn = row.MRN
        total_strains = row["Strain Count"]
        strains_name = row.Strain_x
        total_number_of_mutation_changes = row.mutation_count_x
        total_mutation_changes = row.total_change

        counter=0
        for mut in row.changed_mutations:
            selected_mutation = mut
            counter+=1
            new_mutations.append([mrn,total_strains,strains_name,total_number_of_mutation_changes,total_mutation_changes,str(counter)+"/"+str(total_mutation_changes),selected_mutation])


    df_strain_count_all_mutations = pd.DataFrame(new_mutations)
    df_strain_count_all_mutations.columns=['mrn','total_strains','strain_names','total_number_of_mutation_changes','total_mutation_changes','mutation_counter','selected_mutation']

    df_strain_count_all_mutations["POS"] = [x[1:len(x)-1] for x in df_strain_count_all_mutations.selected_mutation]

    df_strain_count_all_mutations["POS"] = df_strain_count_all_mutations["POS"].astype(int)

    df_strain_count_all_mutations["REF"] = [x[0] for x in df_strain_count_all_mutations.selected_mutation]

    df_strain_count_all_mutations["ALT"] = [x[len(x)-1] for x in df_strain_count_all_mutations.selected_mutation]



    df_strain_count_all_mutations_annotation = pd.merge(df_strain_count_all_mutations,df_snpeff,how="left",on=["POS","REF","ALT"],indicator=True)
    df_strain_count_all_mutations_annotation=df_strain_count_all_mutations_annotation[
    ['mrn','total_strains','strain_names','total_number_of_mutation_changes','total_mutation_changes','mutation_counter','selected_mutation',
           ' Annotation ',' Annotation_Impact ', ' Gene_Name ', ' HGVS.c ',' HGVS.p ']
    ]

    df_strain_count_all_mutations_annotation.to_excel(writer, sheet_name=str(strain_count)+'_strains', index=False)
writer.save()



###################### now analyse this file to generate plots

xls = pd.ExcelFile("Final_snp_analysis_strains_per_patient.xlsx")
df_2_strains = pd.read_excel(xls,'2_strains')
df_3_strains = pd.read_excel(xls,'3_strains')
df_4_strains = pd.read_excel(xls,'4_strains')
df_5_strains = pd.read_excel(xls,'5_strains')
df_6_strains = pd.read_excel(xls,'6_strains')

df_all_strains  = df_2_strains.append(df_3_strains)
df_all_strains  = df_all_strains.append(df_4_strains)
df_all_strains  = df_all_strains.append(df_5_strains)
df_all_strains  = df_all_strains.append(df_6_strains)

df_all_strains.columns  = [x.replace(" ","") for x in df_all_strains.columns]

df_all_strains.to_excel("df_all_strains.xlsx",index=False)

#####
df_all_strains[( (df_all_strains.Gene_Name=="S") & (df_all_strains.Annotation=="missense_variant"))]["HGVS.p"]

df_all_strains[ (df_all_strains.Annotation=="missense_variant")]["Gene_Name"].value_counts()

#################################################


spike = df_all_strains[( (df_all_strains.Gene_Name=="S") & (df_all_strains.Annotation=="missense_variant"))]["HGVS.p"]
df_spike = pd.DataFrame(spike)
df_spike.columns =["AAchange"]
df_spike["AAlocation"] = [ int(re.findall(r'\d+',x)[0]) for x in df_spike["AAchange"].values]
df_spike.sort_values("AAlocation",inplace=True)
df_spike["AAchange2"]=[x.split(".")[1] for x in df_spike.AAchange]
df_spike.drop_duplicates(inplace=True)
df_spike.head()

x= range(1,1273,1)
y =[]
for location in x:
    val=0
    if location in df_spike.AAlocation.values:
        val=1
    y.append(val)

from matplotlib.pyplot import figure
plt.rcParams.update({'font.size': 25})
figure(num=None, figsize=(30, 5), dpi=300, facecolor='w', edgecolor='k')
plt.stem(x,y);

# for location in  df_spike.AAlocation.values:
#     plt.annotate(df_spike[df_spike.AAlocation==location]["AAchange2"].values, (x[location-1], y[location-1]))

plt.yticks([])
plt.ylabel("Mutations observed")
plt.xlabel("Spike protein")
plt.tight_layout()
plt.savefig("temp.png");plt.close()


#################################################

for indx,row in df_mrn_changed_mutation_join.iterrows():

    annotation ={'missense_variant':0,
                'synonymous_variant':0,
                'stop_gained':0,
                'upstream_gene_variant':0,
                'stop_lost&splice_region_variant':0,
                'splice_region_variant&stop_retained_variant':0}

    mut_list = row["changed_mutations"].replace(" ","").replace("[","").replace("]","").replace("'","").split(",")

    for mutation in mut_list:

        if mutation != "":

            pos = int(mutation[1:len(mutation)-1])
            ref = mutation[0]
            alt = mutation[len(mutation)-1]

            new_annotation = df_snpeff[ ((df_snpeff.POS==pos) & (df_snpeff.REF==ref) & (df_snpeff.ALT==alt) )]["Annotation"].values

            annotation[new_annotation[0]] += 1



##### gain or loss over time for patients with 2 strains
df= t2
df = df[df.mutation_diff_1==1]
df2 = df[df["Strain Count"]==2]
df2["res"]=[1 if int(x.split('/')[0])>int(x.split('/')[1]) else 0 for x in df2["mutation_count"]]
df2.res.value_counts()

df2 = df2[df2.res==0] ## for gained
df2 = df2[df2.res==1] ## for loss
df2join = pd.merge(df2,df_mrn_all, on="MRN",how="left")
df2join.groupby("MRN").COLLECTION_DT.diff()
df2join_mini["type"] = df2join.groupby("MRN").COLLECTION_DT.diff().fillna(0)
df2join_mini["mutation_type"] = ["Gain" if x==0 else "Loss" for x in df2join_mini.res]
df2join_mini["interval"] = [x.total_seconds()/(60.0*60.0) for x in df2join_mini.type]
df2join_mini = df2join_mini.drop_duplicates(["MRN"],keep="last")
sns.boxplot(x="mutation_type", y="interval", data=df2join_mini)
sns.swarmplot(x="mutation_type", y="interval", data=df2join_mini, color=".25")
plt.savefig("temp.png");plt.close()

##### for each SNP gained or lost, how many individual times it is gained or lost
df= t2
df = df[df.mutation_diff_1==1]
df2 = df[df["Strain Count"]==2]

df2join = pd.merge(df2,df_mrn_all, on="MRN",how="left")
df2join_mutation = df2join.groupby("MRN").agg( { "mutation_info":'/'.join }).reset_index()

gained ={}
lost ={}
for indx, row in df2join_mutation.iterrows():
    strain_one = row["mutation_info"].split('/')[0]
    strain_two = row["mutation_info"].split('/')[1]

    strain_one = strain_one.replace(" ","").replace("[","").replace("]","").replace("'","").split(",")
    strain_two = strain_two.replace(" ","").replace("[","").replace("]","").replace("'","").split(",")

    if len(strain_one) > len(strain_two): ## loss
        for mutation in strain_one:
            if mutation not in strain_two:
                if mutation not in gained:
                    gained[mutation] = 1
                else:
                    gained[mutation] += 1
    else:
        for mutation in strain_two:
            if mutation not in strain_one:
                if mutation not in lost:
                    lost[mutation] = 1
                else:
                    lost[mutation] += 1

df_gained = pd.DataFrame(list(gained.items()),columns = ['mutation','number_change'])
df_gained["mutation_type"] = "gained"

df_lost = pd.DataFrame(list(lost.items()),columns = ['mutation','number_change'])
df_lost["mutation_type"] = "lost"

df_gained = df_gained.append(df_lost)

df_gained.to_excel("SNP_gained_lost_summary.xlsx",index=False)


######################################
df_curated = pd.read_csv("4_Curated_MCOV_MRN_Strains.csv")
df_curated.columns =['Strain', 'ORDER_ID','MRN','COLLECTION_DT']
df_curated['COLLECTION_DT'] = pd.to_datetime(df_curated['COLLECTION_DT'])

df_mutations_per_strain = pd.read_csv("5_StrainsAndMutations_entire_genome.csv")

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

#####
##### number of changes and strains count

from functools import partial
def checkDiff(mlist,th):
    mutation_count = mlist.split('/')
    first_num = int(mutation_count[0])
    diff = 0
    mut_diff=[]
    for n in mutation_count[1:]:
        mut_diff.append(abs(first_num-int(n)))
    if th==np.max(mut_diff):
        diff = 1
    return diff

for i in range(0,20,1):
    col_name = "mutation_diff_"+str(i)
    t2[col_name] = list(map(partial(checkDiff,th=i),t2.mutation_count.values))


sum_table = t2[["mutation_diff_0","mutation_diff_1","mutation_diff_2","mutation_diff_3","mutation_diff_4",
"mutation_diff_5","mutation_diff_6","mutation_diff_7",
"mutation_diff_8","mutation_diff_9","mutation_diff_10",
"mutation_diff_11","mutation_diff_12","mutation_diff_13",
"mutation_diff_14","mutation_diff_15"
]].sum()

df_sum_table = pd.DataFrame(sum_table)
df_sum_table.columns = ["Total_strains"]
df_sum_table["percent"] =  round((df_sum_table["Total_strains"]/t2.shape[0])*100,2)


#### now look at patients with two strains with one mutation

df_2strain_1mutation = t2[((t2["Strain Count"]==2)& (t2["mutation_diff_1"]==1))]

df_2strain_1mutation = t2[t2["Strain Count"]==2]

df_2strain_1mutation = df_2strain_1mutation[['MRN', 'Strain Count', 'mutation_count', 'Strain']]

df_2strain_1mutation_join = pd.merge(df_2strain_1mutation,df_mrn_all, on="MRN",how="left")

df_2strain_1mutation_join = df_2strain_1mutation_join[['MRN', 'Strain Count', 'mutation_count_x', 'Strain_x', 'ORDER_ID', 'Strain_y', 'mutation_count_y', 'mutation_info']]


mrn_changed_mutation ={}

for mrn in df_2strain_1mutation_join.MRN.unique():

    df_patient = df_2strain_1mutation_join[df_2strain_1mutation_join.MRN==mrn]

    all_mutations = []
    for indx,mut_list in df_patient.mutation_info.items():
        mut_list = ast.literal_eval(mut_list)
        for mut in mut_list:
            if mut not in all_mutations:
                all_mutations.append(mut)

    changed_mutations = []
    for indx,mut_list2 in df_patient.mutation_info.items():
        mut_list2 = ast.literal_eval(mut_list2)
        for mut2 in all_mutations:
            if mut2 not in mut_list2:
                changed_mutations.append(mut2)

    mrn_changed_mutation[mrn] = changed_mutations

df_mrn_changed_mutation = pd.DataFrame(list(mrn_changed_mutation.items()),columns = ['MRN','changed_mutations'])
df_mrn_changed_mutation["total_change"] = [len(x) for x in df_mrn_changed_mutation.changed_mutations]
df_mrn_changed_mutation = df_mrn_changed_mutation[['MRN','total_change','changed_mutations']]
df_mrn_changed_mutation_join = pd.merge(df_2strain_1mutation_join,df_mrn_changed_mutation, on="MRN",how="left")
df_mrn_changed_mutation_join.to_excel("10_mrn_changed_mutation_join.xlsx",index=False)

df_mrn_changed_mutation_join_12 = df_mrn_changed_mutation_join[df_mrn_changed_mutation_join.total_change==1]


df_mrn_changed_mutation_join_12 = df_mrn_changed_mutation_join_12[['MRN', 'Strain Count', 'mutation_count_x', 'Strain_x', 'changed_mutations']]

df_mrn_changed_mutation_join_12["changed_mutations"] = [x[0] for x in df_mrn_changed_mutation_join_12["changed_mutations"] ]

df_mrn_changed_mutation_join_12 = df_mrn_changed_mutation_join_12.drop_duplicates(['MRN','changed_mutations'])

df_mrn_changed_mutation_join_12["POS"] = [x[1:len(x)-1] for x in df_mrn_changed_mutation_join_12.changed_mutations]

df_mrn_changed_mutation_join_12["POS"] = df_mrn_changed_mutation_join_12["POS"].astype(int)

df_mrn_changed_mutation_join_12["REF"] = [x[0] for x in df_mrn_changed_mutation_join_12.changed_mutations]

df_mrn_changed_mutation_join_12["ALT"] = [x[len(x)-1] for x in df_mrn_changed_mutation_join_12.changed_mutations]


VEP_FILE="covid.vcf.snpeff"
header=[]
with open(VEP_FILE) as myfile:
    header = [next(myfile) for x in range(3)]
columns_line=str(header[2].strip().split(':')[1]).split('|')

df_snpeff = pd.read_csv(VEP_FILE,sep='\t',skiprows=5,header=None)
df_snpeff.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
df_snpeff["HGMD_INFO"]= [x.split(";")[0] for x in df_snpeff.INFO]
df_snpeff["VEP_INFO"]= [x.split(";")[1].split(',')[0] for x in df_snpeff.INFO]
df_snpeff[columns_line] = df_snpeff['VEP_INFO'].str.split('|',expand=True)
# df_snpeff = df_snpeff[['CHROM', 'POS', 'REF', 'ALT',' Annotation ']]
# df_snpeff.columns = ['CHROM', 'POS', 'REF', 'ALT','Annotation']

df_mrn_changed_mutation_join_12_annotation_join = pd.merge(df_mrn_changed_mutation_join_12,df_snpeff,how="left",on=["POS","REF","ALT"],indicator=True)

df_mrn_changed_mutation_join_12_annotation_join.to_excel("9_mrn_changed_mutation_join_12_annotation.xlsx",index=False)

###### look at changed mutations

df_mut_1 = t2.copy()
df_mut_1 = df_mut_1[['MRN', 'Strain Count', 'mutation_count', 'Strain']]
df_mut_1_join = pd.merge(df_mut_1,df_mrn_all, on="MRN",how="left")

mrn_changed_mutation ={}

for mrn in df_mut_1_join.MRN.unique():

    df_patient = df_mut_1_join[df_mut_1_join.MRN==mrn]

    all_mutations = []
    for indx,mut_list in df_patient.mutation_info.items():
        mut_list = ast.literal_eval(mut_list)
        for mut in mut_list:
            if mut not in all_mutations:
                all_mutations.append(mut)

    changed_mutations = []
    for indx,mut_list2 in df_patient.mutation_info.items():
        mut_list2 = ast.literal_eval(mut_list2)
        for mut2 in all_mutations:
            if mut2 not in mut_list2:
                changed_mutations.append(mut2)

    mrn_changed_mutation[mrn] = changed_mutations

df_mrn_changed_mutation = pd.DataFrame(list(mrn_changed_mutation.items()),columns = ['MRN','changed_mutations'])
df_mrn_changed_mutation["total_change"] = [len(x) for x in df_mrn_changed_mutation.changed_mutations]
df_mrn_changed_mutation = df_mrn_changed_mutation[['MRN','total_change','changed_mutations']]
df_mrn_changed_mutation_join = pd.merge(df_mut_1,df_mrn_changed_mutation, on="MRN",how="left")
df_mrn_changed_mutation_join.to_excel("8_mrn_changed_mutation_join.xlsx",index=False)

df_mrn_changed_mutation_join = pd.read_excel("8_mrn_changed_mutation_join.xlsx")
####annotate variants
vcf_file = []
for indx,row in df_mrn_changed_mutation_join.iterrows():
    for mutation in row["changed_mutations"]:
        if mutation != "":
            vcf_file.append(["NC_045512.2",mutation[1:len(mutation)-1],".",mutation[0],mutation[len(mutation)-1],"100.0","PASS","INFO"])
df_vcf_file = pd.DataFrame(vcf_file)
df_vcf_file[1] = df_vcf_file[1].astype(int)
df_vcf_file = df_vcf_file.drop_duplicates()
df_vcf_file = df_vcf_file.sort_values(1)
df_vcf_file.to_csv("covid.vcf",sep='\t',header=None,index=False)


VEP_FILE="covid.vcf.snpeff"
header=[]
with open(VEP_FILE) as myfile:
    header = [next(myfile) for x in range(3)]
columns_line=str(header[2].strip().split(':')[1]).split('|')

df_snpeff = pd.read_csv(VEP_FILE,sep='\t',skiprows=5,header=None)
df_snpeff.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
df_snpeff["HGMD_INFO"]= [x.split(";")[0] for x in df_snpeff.INFO]
df_snpeff["VEP_INFO"]= [x.split(";")[1].split(',')[0] for x in df_snpeff.INFO]
df_snpeff[columns_line] = df_snpeff['VEP_INFO'].str.split('|',expand=True)
df_snpeff = df_snpeff[['CHROM', 'POS', 'REF', 'ALT',' Annotation ']]
df_snpeff.columns = ['CHROM', 'POS', 'REF', 'ALT','Annotation']


for indx,row in df_mrn_changed_mutation_join.iterrows():

    annotation ={'missense_variant':0,
                'synonymous_variant':0,
                'stop_gained':0,
                'upstream_gene_variant':0,
                'stop_lost&splice_region_variant':0,
                'splice_region_variant&stop_retained_variant':0}

    mut_list = row["changed_mutations"].replace(" ","").replace("[","").replace("]","").replace("'","").split(",")

    for mutation in mut_list:

        if mutation != "":

            pos = int(mutation[1:len(mutation)-1])
            ref = mutation[0]
            alt = mutation[len(mutation)-1]

            new_annotation = df_snpeff[ ((df_snpeff.POS==pos) & (df_snpeff.REF==ref) & (df_snpeff.ALT==alt) )]["Annotation"].values

            annotation[new_annotation[0]] += 1



##### gain or loss over time for patients with 2 strains
df= t2
df = df[df.mutation_diff_1==1]
df2 = df[df["Strain Count"]==2]
df2["res"]=[1 if int(x.split('/')[0])>int(x.split('/')[1]) else 0 for x in df2["mutation_count"]]
df2.res.value_counts()

df2 = df2[df2.res==0] ## for gained
df2 = df2[df2.res==1] ## for loss
df2join = pd.merge(df2,df_mrn_all, on="MRN",how="left")
df2join.groupby("MRN").COLLECTION_DT.diff()
df2join_mini["type"] = df2join.groupby("MRN").COLLECTION_DT.diff().fillna(0)
df2join_mini["mutation_type"] = ["Gain" if x==0 else "Loss" for x in df2join_mini.res]
df2join_mini["interval"] = [x.total_seconds()/(60.0*60.0) for x in df2join_mini.type]
df2join_mini = df2join_mini.drop_duplicates(["MRN"],keep="last")
sns.boxplot(x="mutation_type", y="interval", data=df2join_mini)
sns.swarmplot(x="mutation_type", y="interval", data=df2join_mini, color=".25")
plt.savefig("temp.png");plt.close()

##### for each SNP gained or lost, how many individual times it is gained or lost
df= t2
df = df[df.mutation_diff_1==1]
df2 = df[df["Strain Count"]==2]

df2join = pd.merge(df2,df_mrn_all, on="MRN",how="left")
df2join_mutation = df2join.groupby("MRN").agg( { "mutation_info":'/'.join }).reset_index()

gained ={}
lost ={}
for indx, row in df2join_mutation.iterrows():
    strain_one = row["mutation_info"].split('/')[0]
    strain_two = row["mutation_info"].split('/')[1]

    strain_one = strain_one.replace(" ","").replace("[","").replace("]","").replace("'","").split(",")
    strain_two = strain_two.replace(" ","").replace("[","").replace("]","").replace("'","").split(",")

    if len(strain_one) > len(strain_two): ## loss
        for mutation in strain_one:
            if mutation not in strain_two:
                if mutation not in gained:
                    gained[mutation] = 1
                else:
                    gained[mutation] += 1
    else:
        for mutation in strain_two:
            if mutation not in strain_one:
                if mutation not in lost:
                    lost[mutation] = 1
                else:
                    lost[mutation] += 1

df_gained = pd.DataFrame(list(gained.items()),columns = ['mutation','number_change'])
df_gained["mutation_type"] = "gained"

df_lost = pd.DataFrame(list(lost.items()),columns = ['mutation','number_change'])
df_lost["mutation_type"] = "lost"

df_gained = df_gained.append(df_lost)

df_gained.to_excel("SNP_gained_lost_summary.xlsx",index=False)
