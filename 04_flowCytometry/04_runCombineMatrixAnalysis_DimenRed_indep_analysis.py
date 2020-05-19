import pandas as pd
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import pylab
import pandas as pd
import seaborn as sns
import scipy.stats as stats
from sklearn import decomposition
from sklearn.manifold import TSNE


def concatRows(row):
    return str(row[0])+str(row[1])+str(row[2])+str(row[3])+str(row[4])+str(row[5])+str(row[6])


df_sample = pd.read_csv("/net/fs01.cluster.com/home/hhadmin/flowCyto/indepen_analysis/results/t_cells_result/FLW195056_tcells_matrix.csv")



df_sample["combine"] = df_sample[['PE.ASSC.APOS_NEG', 'PerCP.Cy5.5.AAPC.H7.APOS_NEG','PerCP.Cy5.5.AAPC.H7.ANEG_POS', 'APC.ASSC.APOS_POS',
'APC.ASSC.APOS_NEG', 'PE.AAPC.APOS_POS', 'PE.AAPC.ANEG_POS']].apply(concatRows,axis=1)
# selected = ['1010000','1100000','1100110','0000101','1000000','0000001','0001001','1000110','0010000','0100000']
selected = ['1010000','1100000','1100110','0000101','0001001','1000110']
df_sample = df_sample[df_sample["combine"].isin(selected)]
df_sample["combine"] = ["$"+x for x in df_sample["combine"].values ]

####PCA analysis
### combination analysis

X1= df_sample.iloc[:,7:15].values
pca = decomposition.PCA(n_components=4)
pc = pca.fit_transform(X1)

pc_df = pd.DataFrame(data = pc ,columns = ['PC1', 'PC2','PC3','PC4'])
pc_df['Cluster'] = df_sample["combine"].values
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
  hue='Cluster', # color by cluster
  legend=True,
  scatter_kws={"s": 4}) # specify the point size
plt.savefig("pca_2.png");plt.close()

############ tsne

X1= df_sample.iloc[:,7:15].values
tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
tsne_results = tsne.fit_transform(X1)
df_subset = pd.DataFrame(data = tsne_results ,columns = ['tsne-2d-one', 'tsne-2d-two'])
df_subset['Cluster'] = df_sample["combine"].values
plt.figure(figsize=(16,10))
sns.scatterplot(
    x="tsne-2d-one", y="tsne-2d-two",
    hue="Cluster",
    palette=sns.color_palette("hls", 5),
    data=df_subset,
    legend="full",
)
plt.savefig("tsne_1.png");plt.close()

##### heatmap and boxplot analysis

df= df_sample.iloc[:,7:15]
sns.heatmap(df,cmap="YlGnBu")
plt.savefig("heatmap_1.png");plt.close()

sns.boxplot(x="variable", y="value", data=pd.melt(df))
plt.savefig("boxplot_1.png");plt.close()



##################### b cells


def concatRowsB(row):
    return str(row[0])+str(row[1])+str(row[2])+str(row[3])+str(row[4])+str(row[5])+str(row[6])+str(row[7])


df_sample = pd.read_csv("/net/fs01.cluster.com/home/hhadmin/flowCyto/indepen_analysis/results/b_cells_result/FLW195056_bcells_matrix.csv")



df_sample["combine"] = df_sample[
        ['V450.ASSC.APOS_NEGB', 'PE.AFITC.APOS_NEGB-L',
       'PE.AFITC.ANEG_POSB-K', 'V450.APerCP.Cy5.5.APOS_POSB-CD5',
       'APC.ASSC.APOS_NEGCD10', 'PE.AFITC.APOS_NEGCD10-L',
       'PE.AFITC.ANEG_POSCD10-K', 'APC.H7.ASSC.APOS_POSPLASMA']].apply(concatRowsB,axis=1)
selected = ['10100000','00001000','00000001','11010000','10001000','00010000','00001001','10001001']
df_sample = df_sample[df_sample["combine"].isin(selected)]
df_sample["combine"] = ["$"+x for x in df_sample["combine"].values ]

####PCA analysis
### combination analysis

X1= df_sample.iloc[:,7:15].values
pca = decomposition.PCA(n_components=4)
pc = pca.fit_transform(X1)

pc_df = pd.DataFrame(data = pc ,columns = ['PC1', 'PC2','PC3','PC4'])
pc_df['Cluster'] = df_sample["combine"].values
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
  hue='Cluster', # color by cluster
  legend=True,
  scatter_kws={"s": 4}) # specify the point size
plt.savefig("pca_2.png");plt.close()

############ tsne

X1= df_sample.iloc[:,7:15].values
tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
tsne_results = tsne.fit_transform(X1)
df_subset = pd.DataFrame(data = tsne_results ,columns = ['tsne-2d-one', 'tsne-2d-two'])
df_subset['Cluster'] = df_sample["combine"].values
plt.figure(figsize=(16,10))
sns.scatterplot(
    x="tsne-2d-one", y="tsne-2d-two",
    hue="Cluster",
    palette=sns.color_palette("hls", 8),
    data=df_subset,
    legend="full",
)
plt.savefig("tsne_1.png");plt.close()

##### heatmap and boxplot analysis

df= df_sample.iloc[:,7:15]
sns.heatmap(df,cmap="YlGnBu")
plt.savefig("heatmap_1.png");plt.close()

sns.boxplot(x="variable", y="value", data=pd.melt(df))
plt.savefig("boxplot_1.png");plt.close()
