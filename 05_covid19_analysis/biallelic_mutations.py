import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import seaborn as sns
from Bio import AlignIO, SeqIO, Entrez
import numpy as np
import os
import re
import subprocess
import ast
import vcf_analysis


############### biallelic_mutations
vcfs_directory="data/7_15/all_vcfs/"
df_vcf = vcf_analysis.analyzeVCF(vcfs_directory)
df_vcf.to_csv("reference/all_strains_vcf_record.csv",index=False)

df=pd.read_excel("data/biallelic_mutations/biallelicmutations_nsp12.xlsx")


filter =[]
for indx,row in df.iterrows():
    df_mini = df_vcf[df_vcf["POS"]==row.NTA_Genomic_Locus]

    vdict = ast.literal_eval(row.NTA_Variants_Info)

    strains =[]
    for items in vdict:
        for strain in vdict[items][1]:
            strains.append(strain)

    for strain in strains:
        if df_mini[df_mini.Strain==strain].shape[0] ==1:
            strain_info = df_mini[df_mini.Strain==strain].values[0]
            filter.append([row.NTA_Genomic_Locus, row.NTA_Reference.upper(),strain_info[0],strain_info[3],strain_info[4]])

df_bialleleic_snps = pd.DataFrame(filter, columns=["Position","Reference","Strain","Variant","Quality"])
df_bialleleic_snps.to_excel("df_bialleleic_snps_nsp12.xlsx",index=False)

dfjoin = pd.merge(df,df_bialleleic_snps, how="left",right_on="Position",left_on="NTA_Genomic_Locus",indicator=True)
dfjoin.to_excel("bialleleic_nta_vcf_comparison_nsp12.xlsx",index=False)
################


# df=pd.read_excel("data/biallelic_mutations/COVID_SNP_MAP_S_NTA_7_19.xlsx")
df=pd.read_excel("data/biallelic_mutations/COVID_SNP_MAP_nsp12_NTA_7_19.xlsx")

df = df[['NTA_Genomic_Locus', 'NTA_Reference','NTA_Variants','Annotation','HGVS.p.short']]
df["NTA_Reference"] = [x.upper() for x in df["NTA_Reference"]]

nvariant=[]
for indx,row in df.iterrows():
    vlist = ast.literal_eval(row.NTA_Variants)
    vlen= len(vlist)
    for v in vlist:
        nvariant.append([row.NTA_Genomic_Locus,row.NTA_Reference,v.upper(),row.Annotation,row["HGVS.p.short"],vlen])

df2 = pd.DataFrame(nvariant, columns=['NTA_Genomic_Locus', 'NTA_Reference','NTA_Variants','Annotation','HGVS.p.short','len'])



###########################################
df_db= pd.read_excel("data/biallelic_mutations/Database_Variation_Annotation_8_17.xlsx")
df_db = df_db[['Genome position', 'Base change:Virus number']]


df=pd.read_excel("data/biallelic_mutations/COVID_SNP_MAP_S_NTA_table_8_11.xlsx")
# df=pd.read_excel("data/biallelic_mutations/COVID_SNP_MAP_nsp12_NTA_table_8_11.xlsx")

dfjoin = pd.merge(df,df_db, how="left",right_on="Genome position",left_on="Genomic locus",indicator=True)
# dfjoin["status"] = [x.upper()in(str(y))  for x,y  in zip(dfjoin.NTA_Variants,dfjoin["Base change:Virus number"])]

position_dict = dfjoin["Genomic locus"].value_counts().to_dict()
dfjoin["biallelic"] = [position_dict[x] for x in dfjoin["Genomic locus"]]

# dfjoin.to_excel("data/biallelic_mutations/bialleleic_nta_cndb_comparison_all_s_protein.xlsx",index=False)
dfjoin.to_excel("data/biallelic_mutations/bialleleic_nta_cndb_comparison_all_nsp12_protein.xlsx",index=False)



# ### get snp figure

df=pd.read_excel("data/biallelic_mutations/COVID_SNP_MAP_S_NTA_table_8_11.xlsx")

total_samples = 5085.0
aa_length=1274
aa_cordinates=range(1,aa_length,1)
mut_dict = {x:0 for x in aa_cordinates}
for indx,row in df.iterrows():
    aa_loc = row["AA Change"]
    if aa_loc != "Syn":
        position =  int(re.sub('[A-Z]','',aa_loc))
        # mut_proportion = float(row['Total Confirmed (5085)'])/total_samples
        mut_proportion = np.log10(row['Total Confirmed (5085)'])
        mut_proportion = row['Total Confirmed (5085)']
        # print(row['Total Confirmed (5085)'],np.log(row['Total Confirmed (5085)']),np.log10(row['Total Confirmed (5085)']))
        if mut_dict[position] != 0:
            mut_dict[position] += mut_proportion
        else:
            mut_dict[position] = mut_proportion



from matplotlib.pyplot import figure
figure(num=None, figsize=(10, 4), dpi=300, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 15})

sns.lineplot(aa_cordinates, y=[x for x in mut_dict.values()],color='black') ##remove first sampleID to plot

plt.title("Isolates (N=5085)")
plt.ylabel("Isolates Containing SNP")
plt.xlabel("S protein AA codon")
plt.tight_layout()
plt.savefig("COVID_SNP_MAP_S_protein_v2_raw.png");
plt.close()
