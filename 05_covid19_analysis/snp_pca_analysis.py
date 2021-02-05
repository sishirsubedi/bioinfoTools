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
from sklearn import decomposition


df_db = pd.read_csv("1_StrainsAndMutations_entire_genome.csv")
df_db = df_db[df_db.mutation_count!=0]


mutations = []
for indx,row in df_db.iterrows():
    mut_list = ast.literal_eval(row.mutation_info)
    for m in mut_list:
        if m not in mutations:
            mutations.append(m)

##if we want to order muations
df_mut = pd.DataFrame(mutations)
df_mut.columns=["mutations"]
df_mut["location"] = [x[1:len(x)-1] for x in df_mut.mutations]
df_mut["location"] = df_mut["location"].astype(int)
df_mut.sort_values("location",inplace=True)


strain_mutations =[]
for indx,row in df_db.iterrows():

    strain=row.Strain
    mut_list = ast.literal_eval(row.mutation_info)

    new_row =[]
    new_row.append(strain)

    for m in mutations:
        if m in mut_list:
            new_row.append(1)
        else:
            new_row.append(0)

    strain_mutations.append(new_row)

df_strain_mutations = pd.DataFrame(strain_mutations)
new_cols = []
new_cols.append("strain")
for m in mutations: new_cols.append(m)
df_strain_mutations.columns = new_cols

df_strain_mutations.to_excel("out_1_df_strain_mutations_1_0_matrix.xlsx",index=False)



###############################
### pca analysis #############
#############################


X1= df_strain_mutations.iloc[:,1:].values
pca = decomposition.PCA(n_components=4)
pc = pca.fit_transform(X1)

pc_df = pd.DataFrame(data = pc ,columns = ['PC1', 'PC2','PC3','PC4'])
pc_df['D614G'] = ["Absent" if x==0 else "Present" for x in df_strain_mutations.A23403G]
# pc_df['D614G'] = ["Absent" if x==0 else "Present" for x in df_strain_mutations.C14408T]
# pc_df['D614G'] = ["Absent" if x==0 else "Present" for x in df_strain_mutations.C3037T]


pc_df.head()

pca.explained_variance_ratio_
df = pd.DataFrame({'var':pca.explained_variance_ratio_,
             'PC':['PC1','PC2','PC3','PC4']})
sns.barplot(x='PC',y="var",
           data=df, color="c");
plt.savefig("pca_1.png");plt.close()

sns.lmplot( x="PC1", y="PC2",
  data=pc_df,
  fit_reg=False,
  hue='D614G', # color by cluster
  legend=True,
  scatter_kws={"s": 4}) # specify the point size
plt.savefig("pca_2.png");plt.close()


sns.heatmap(df_strain_mutations.iloc[:,1:].corr(),vmin=-1,
            cmap='coolwarm',
            annot=True);
plt.savefig("heatmap.png");plt.close()
