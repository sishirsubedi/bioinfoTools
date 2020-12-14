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

def concatRows(cell,row):
    if cell=="T":
        return str(row[0])+str(row[1])+str(row[2])+str(row[3])+str(row[4])+str(row[5])+str(row[6])
    elif cell =="B":
        return str(row[0])+str(row[1])+str(row[2])+str(row[3])+str(row[4])+str(row[5])+str(row[6])+str(row[7])

def processTcells(result_file,result_dir):

    df = pd.read_csv(result_file,header=None)
    df.columns =['AccessionNumber','Marker','cell_count','cell_proportion_auto','flow_index']
    df["Marker"] = [x.replace("[","").replace("]","").replace(" ","") for x in df.Marker.values]
    df['AccessionNumber'] = [ x.replace(" ","") for x in df.AccessionNumber.values ]


    samples = ['FLW195056'] ## all samples of passed only

    for sample in samples:
        print("processing...."+sample)
        df_sample = pd.read_csv( result_dir+" "+sample+" _trans_mod.csv",header=None)
        df_sample.columns = ["Time","FSC.A","FSC.H","FSC.W","SSC.A","SSC.H","SSC.W","FITC.A","PE.A","PerCP.Cy5.5.A","PE.Cy7.A","APC.A","APC.H7.A","V450.A","V500.A"]
        df_result = df.ix[df.AccessionNumber==sample,]

        for indx,row in df_result.iterrows():
            df_sample[row["Marker"]] = 0
            for findex in row["flow_index"].split(";"):
                df_sample.ix[int(findex.replace(" ",""))-1,row["Marker"]] = 1

        df_sample["combine"] = df_sample[['PE.ASSC.APOS_NEG', 'PerCP.Cy5.5.AAPC.H7.APOS_NEG','PerCP.Cy5.5.AAPC.H7.ANEG_POS', 'APC.ASSC.APOS_POS',
        'APC.ASSC.APOS_NEG', 'PE.AAPC.APOS_POS', 'PE.AAPC.ANEG_POS']].apply(lambda x : concatRows("T",x),axis=1)
        selected = ['1010000','1100000','1100110','0000101']
        df_sample = df_sample[df_sample["combine"].isin(selected)]

        markers ={'1010000' :'CD4',
        '1100000':'CD8',
        '1100110': 'NK like T LGL',
        '0000101': 'NK'}

        df_sample["combine2"] = [ markers[x] for x in df_sample["combine"].values]
        # df_sample["combine2"] = [ x for x in df_sample["combine"].values]

        df_sample.sort_values(by="combine2",inplace=True)

        X1= df_sample.iloc[:,7:15].values
        pca = decomposition.PCA(n_components=4)
        pc = pca.fit_transform(X1)

        pc_df = pd.DataFrame(data = pc ,columns = ['PC1', 'PC2','PC3','PC4'])
        pc_df['Cluster'] = df_sample["combine2"].values
        plt.figure(figsize=(16,10))
        sns.lmplot( x="PC1", y="PC2",
        data=pc_df,
        fit_reg=False,
        hue='Cluster', # color by cluster
        legend=False,
        scatter_kws={"s": 4}) # specify the point size
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel("PCA1",fontsize=20)
        plt.ylabel("PCA2",fontsize=20)
        plt.legend(prop={'size': 20},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+sample+"_PCA_1.png",bbox_inches = "tight",dpi=300)
        plt.close()

def processBcells(result_file,result_dir):

    df = pd.read_csv(result_file,header=None)
    df.columns =['AccessionNumber','Marker','cell_count','cell_proportion_auto','flow_index']
    df["Marker"] = [x.replace("[","").replace("]","").replace(" ","") for x in df.Marker.values]
    df['AccessionNumber'] = [ x.replace(" ","") for x in df.AccessionNumber.values ]


    samples = ['FLW195056'] ## pass all samples or passed only

    for sample in samples:
        print("processing...."+sample)
        df_sample = pd.read_csv( result_dir+" "+sample+" _trans_mod.csv",header=None)
        df_sample.columns = ["Time","FSC.A","FSC.H","FSC.W","SSC.A","SSC.H","SSC.W","FITC.A","PE.A","PerCP.Cy5.5.A","PE.Cy7.A","APC.A","APC.H7.A","V450.A","V500.A"]
        df_result = df.ix[df.AccessionNumber==sample,]

        for indx,row in df_result.iterrows():
            df_sample[row["Marker"]] = 0
            for findex in row["flow_index"].split(";"):
                if findex == '' or findex == ' ':
                    continue
                else:
                    df_sample.ix[int(findex.replace(" ",""))-1,row["Marker"]] = 1

        df_sample["combine"] = df_sample[['V450.ASSC.APOS_NEGB', 'PE.AFITC.APOS_NEGB-L',
        'PE.AFITC.ANEG_POSB-K', 'V450.APerCP.Cy5.5.APOS_POSB-CD5',
        'APC.ASSC.APOS_NEGCD10', 'PE.AFITC.APOS_NEGCD10-L',
        'PE.AFITC.ANEG_POSCD10-K', 'APC.H7.ASSC.APOS_POSPLASMA']].apply(lambda x : concatRows("B",x),axis=1)


        selected = ['10010000', '10100000', '11000000']

        df_sample = df_sample[df_sample["combine"].isin(selected)]

        markers = {
            "10010000"    : "CD5+ B cell",
            "10100000"    : "B cell (Kappa)",
            "11000000"    : "B cell (Lambda)",
            }

        df_sample["combine2"] = [ markers[x] for x in df_sample["combine"].values]
        # df_sample["combine2"] = [ x for x in df_sample["combine"].values]

        df_sample.sort_values(by="combine2",inplace=True)

        X1= df_sample.iloc[:,7:15].values
        pca = decomposition.PCA(n_components=4)
        pc = pca.fit_transform(X1)

        pc_df = pd.DataFrame(data = pc ,columns = ['PC1', 'PC2','PC3','PC4'])
        pc_df['Cluster'] = df_sample["combine2"].values
        plt.figure(figsize=(16,10))
        sns.lmplot( x="PC1", y="PC2",
        data=pc_df,
        fit_reg=False,
        hue='Cluster', # color by cluster
        legend=False,
        scatter_kws={"s": 50}) # specify the point size
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel("PCA1",fontsize=20)
        plt.ylabel("PCA2",fontsize=20)
        plt.legend(prop={'size': 20},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+sample+"_PCA_1.png",bbox_inches = "tight",dpi=300)
        plt.close()

cell_type = sys.argv[1]
result_file = sys.argv[2]
result_dir = sys.argv[3]

if cell_type == "Tcells":
    processTcells(result_file,result_dir)
elif cell_type =="Bcells":
    processBcells(result_file,result_dir)
