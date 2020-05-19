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
import numpy as np
from sklearn import decomposition
from sklearn.manifold import TSNE

def concatRows(row):
    return str(row[0])+str(row[1])+str(row[2])+str(row[3])+str(row[4])+str(row[5])+str(row[6])

result_file = "/net/fs01.cluster.com/home/hhadmin/flowCyto/data/T_CELLS_RESULT_OUTPUT_TCells_completed/flow_results/ T_CELLS_result.csv"
result_dir = "/net/fs01.cluster.com/home/hhadmin/flowCyto/data/T_CELLS_RESULT_OUTPUT_TCells_completed/"

df = pd.read_csv(result_file,header=None)
df.columns =['AccessionNumber','Marker','cell_count','cell_proportion_auto','flow_index']
df["Marker"] = [x.replace("[","").replace("]","").replace(" ","") for x in df.Marker.values]
df['AccessionNumber'] = [ x.replace(" ","") for x in df.AccessionNumber.values ]


# samples = ['FLW19110',
# 'FLW1919',
# 'FLW1933',
# 'FLW1935',
# 'FLW194503',
# 'FLW194520',
# 'FLW194551',
# 'FLW194555',
# 'FLW194557',
# 'FLW194591',
# 'FLW194597',
# 'FLW194735',
# 'FLW194736',
# 'FLW194815',
# 'FLW194853',
# 'FLW194858',
# 'FLW194876',
# 'FLW194910',
# 'FLW194911',
# 'FLW194971',
# 'FLW194977',
# 'FLW194991',
# 'FLW195004',
# 'FLW195028',
# 'FLW195041',
# 'FLW195043',
# 'FLW195056',
# 'FLW195091']


samples = [
'FLW1935',
'FLW194551',
'FLW194591',
'FLW195056']

for sample in samples:
    print("processing...."+sample)
    df_sample = pd.read_csv( result_dir+"flow_results/"+" "+sample+" _trans_mod.csv",header=None)
    df_sample.columns = ["Time","FSC.A","FSC.H","FSC.W","SSC.A","SSC.H","SSC.W","FITC.A","PE.A","PerCP.Cy5.5.A","PE.Cy7.A","APC.A","APC.H7.A","V450.A","V500.A"]
    df_result = df.ix[df.AccessionNumber==sample,]

    for indx,row in df_result.iterrows():
        df_sample[row["Marker"]] = 0
        for findex in row["flow_index"].split(";"):
            df_sample.ix[int(findex.replace(" ",""))-1,row["Marker"]] = 1

    df_sample["combine"] = df_sample[['PE.ASSC.APOS_NEG', 'PerCP.Cy5.5.AAPC.H7.APOS_NEG','PerCP.Cy5.5.AAPC.H7.ANEG_POS', 'APC.ASSC.APOS_POS',
    'APC.ASSC.APOS_NEG', 'PE.AAPC.APOS_POS', 'PE.AAPC.ANEG_POS']].apply(concatRows,axis=1)
    selected = ['1010000','1100000','1100110','0000101']
    df_sample = df_sample[df_sample["combine"].isin(selected)]

    markers ={'1010000' :'CD4',
    '1100000':'CD8',
    '1100110': 'NK like T LGL',
    '0000101': 'NK'}

    df_sample["combine2"] = [ markers[x] for x in df_sample["combine"].values]
    df_sample.sort_values(by="combine2",inplace=True)
    
    ############ tsne

    X1= df_sample.iloc[:,7:15].values
    tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300,random_state=np.random.RandomState(0))
    tsne_results = tsne.fit_transform(X1)
    df_subset = pd.DataFrame(data = tsne_results ,columns = ['tsne-2d-one', 'tsne-2d-two'])
    df_subset['Cluster'] = df_sample["combine2"].values
    plt.figure(figsize=(16,10))
    sns.scatterplot(
        x="tsne-2d-one", y="tsne-2d-two",
        hue="Cluster",
        palette=['green','red','dodgerblue','orange'],
        data=df_subset,
        legend="full",
    )
    plt.legend(prop={'size': 12},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
    plt.savefig(result_dir+"tsne_plots/"+sample+"_tsne_1.png",bbox_inches = "tight",dpi=300);plt.close()
