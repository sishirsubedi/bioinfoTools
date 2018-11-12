hre
# 1. convert database such as clinvar or cosmic in vcf file formate to appropriate format for running bedtools
# desired format is:
#                       chr      start    end     id              gene
#                       chr1    69224   69225   COSM3677745     OR4F5

awk -v OFS=$'\t' '{ $1="chr"$1;print}' db_02_cosmic_modified.vcf |
awk '{print $1"\t"$2"\t" $2+1 "\t" $0 }' |
grep 'GENE' |
cut -f -3,6,11 |
awk -F';' '{sub("GENE=", "", $1); print $1  }'' |
awk '{split($5,a,"_"); print $1 "\t" $2 "\t" $3  "\t" $4 "\t"  a[1] }' > db_02_cosmic_modified2.vcf


#run bedtools and get required information i.e. gene or postion
/opt/bedtools2/bin/bedtools intersect -a sample_mod.bed -b db_01_clinvar_modified2.vcf -wa -wb | awk '{print $8}' > out_clinvar_v7_intersect_genes.csv


#################################################
#### use python pandas to get genes
#################################################

import pandas as pd

df_db = pd.read_csv("db_01_clinvar_modified2.vcf", sep='\t',header=None)
print(df_db.shape)
df_db.columns = ['chromosome', 'position1', 'position2', 'clinvar_id', 'gene' , 'exon']
filter_indx = df_db[df_db['chromosome'].str.contains("MT|NW_00")].index
df_db.drop(filter_indx, inplace=True)
print(df_db.shape)
df_db['exon_count'] = df_db.groupby('gene')['gene'].transform('count')
df_db2 = df_db.drop_duplicates('gene').reset_index()
df_db2 = df_db2.drop(['index','exon'], 1)
print(df_db2.head())
print(df_db2.shape)

df_db = pd.read_csv("db_03_1_UCSC_genes_coordinates_modified2.vcf", sep='\t',header=None)
print(df_db.shape)
df_db.columns = ['chromosome', 'position1', 'position2', 'refseq_id', 'gene']
df_db['exon_count'] = df_db.groupby('gene')['gene'].transform('count')
df_db2 = df_db.drop_duplicates('gene').reset_index()
df_db2 = df_db2.drop(['index','position1', 'position2','refseq_id','exon_count'], 1)
print(df_db2.head())
print(df_db2.shape)

df_db = pd.read_csv("db_03_1_UCSC_genes_coordinates_modified3_li_status.vcf", sep='\t',header=None)
print(df_db.shape)
df_db.columns = ['chromosome', 'position1', 'position2', 'refseq_id', 'gene' , 'status']
df_db['exon_count'] = df_db.groupby('gene')['gene'].transform('count')
df_db2 = df_db.drop_duplicates('gene').reset_index()
df_db2 = df_db2.drop(['index','position1', 'position2','refseq_id','exon_count'], 1)
print(df_db2.head())
print(df_db2.shape)

df_v7 = pd.read_csv("/home/hhadmin/agilent/v7_cre_comparision/in_design_files/v7/out_ucsc_v7_intersect_genes.csv",header=None)
print(df_v7.head())
v7_genes = df_v7.iloc[:,0].unique()
print(len(v7_genes))

df_cre = pd.read_csv("/home/hhadmin/agilent/v7_cre_comparision/in_design_files/cre/out_ucsc_cre_intersect_genes.csv",header=None)
print(df_cre.head())
cre_genes = df_cre.iloc[:,0].unique()
print(len(cre_genes))

df_twistbio = pd.read_csv("/home/hhadmin/agilent/twistbio/out_ucsc_twistbio_intersect_genes.csv",header=None)
print(df_twistbio.head())
twistbio_genes = df_twistbio.iloc[:,0].unique()
print(len(twistbio_genes))

#v7_unique = [gene for gene in v7_genes if gene not in cre_genes]
df_db2['twistbio'] = df_db2.apply(lambda row: 1 if row['gene'] in twistbio_genes else 0,axis=1)
df_db2['v7'] = df_db2.apply(lambda row: 1 if row['gene'] in v7_genes else 0,axis=1)
df_db2['cre'] = df_db2.apply(lambda row: 1 if row['gene'] in cre_genes else 0, axis=1)

print(df_db2.head())

df_db2.to_csv("out_ucsc_comparison_cre_v7_twistbio_design.csv",index=False)

###################################
### li-list get coordinates from ucsc
df_db = pd.read_csv("db_03_1_UCSC_genes_coordinates_modified2.vcf", sep='\t',header=None)
print(df_db.shape)
df_db.columns = ['chromosome', 'position1', 'position2', 'refseq_id', 'gene']
filter_indx = df_db[df_db['chromosome'].str.contains("_")].index
df_db.drop(filter_indx, inplace=True)

df_genelist = pd.read_csv("db_03_1_Li_Genelist.txt")
df_genelist.head()
df_genelist.drop_duplicates(["gene","type"],inplace=True)
genes = df_genelist.iloc[:,0].unique()

status=[]
for indx,row in df_db.iterrows():
    if row['gene'] in genes:
        status.append(df_genelist[df_genelist['gene']==row['gene']]['type'].values[0])
    else:
        status.append('unknown')
df_db['status'] = status

df_db.to_csv('db_03_1_UCSC_genes_coordinates_modified3_li_status.vcf',index=False,sep='\t',header=None)

###################################
import pandas as pd

healthy = pd.read_csv("healthy_genelist.csv",header=None)
healthy.head()
tumor=pd.read_csv("tumor_genelist.csv",header=None)
tumor.shape

gene_list=[]

for indx,row in healthy.iterrows():
  for v in row:
    if v =='NaN' or not v:
      continue
    else:
      gene_list.append([v,"healthy"])

for indx,row in tumor.iterrows():
  for v in row:
    if v =='NaN' or not v:
      continue
    else:
      gene_list.append([v,"tumor"])
df_genelist = pd.DataFrame(gene_list)
df_genelist.columns=['gene','type']
df_genelist = df_genelist.dropna()
df_genelist.to_csv("genelist_all.csv",index=False)

##################################################


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

venn2(subsets = (56,4, 49), set_labels = ('v7', 'cre'))
plt.savefig('clinvar_bam')
plt.close()

venn2(subsets = (126,120,144), set_labels = ('v7', 'cre'))
plt.savefig('clinvar_design')
plt.close()



import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

 flags = pd.read_csv("samflags",header=None)

 flags = pd.read_csv("mapquals",header=None)

 flags.columns = ['flag']
 flags['count'] = flags.groupby('flag')['flag'].transform('count')
 flag2 = flags.drop_duplicates('flag').reset_index()
 flag2 = flag2.drop(['index'], 1)
 flag2 = flag2.sort_values('flag')
 print(flag2.head())
 print(flag2.shape)

plt.rcParams["figure.figsize"] = (10,10)

flag2.plot.bar(x='flag',y='count',ylim=(0,1000));plt.savefig('test2');plt.close()
