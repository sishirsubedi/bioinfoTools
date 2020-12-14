import pandas as pd
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import pandas as pd
import seaborn as sns
import scipy.stats as stats


def concatRows(row):
    return str(row[0])+str(row[1])+str(row[2])+str(row[3])+str(row[4])+str(row[5])+str(row[6])

def concatRowsBcells(row):
    return str(row[0])+str(row[1])+str(row[2])+str(row[3])+str(row[4])+str(row[5])+str(row[6])+str(row[7])

def plotTcells(result_file,result_dir):

    df = pd.read_csv(result_file,header=None)
    df.columns =['AccessionNumber','Marker','cell_count','cell_proportion_auto','flow_index']
    df["Marker"] = [x.replace("[","").replace("]","").replace(" ","") for x in df.Marker.values]
    df['AccessionNumber'] = [ x.replace(" ","") for x in df.AccessionNumber.values ]

    ####################################### matrix analysis


    cell_marker = {
           'PE.ASSC.APOS_NEG': 'CD3+',
           'PerCP.Cy5.5.AAPC.H7.ANEG_POS':  'CD4+CD8-',
           'PerCP.Cy5.5.AAPC.H7.APOS_NEG': 'CD8+CD4-' ,
           'APC.ASSC.APOS_POS' : 'CD56+SSC+',
           'APC.ASSC.APOS_NEG' : 'CD56+SSC-',
           'PE.AAPC.APOS_POS' :  'CD3+SSC-',
           'PE.AAPC.ANEG_POS' : 'CD3+SSC-'
           }

    for sample in df.AccessionNumber.unique():

        print("processing..."+sample)

        df_sample = pd.read_csv( result_dir+" "+sample+" _trans_mod.csv",header=None)
        df_sample.columns = ["Time","FSC.A","FSC.H","FSC.W","SSC.A","SSC.H","SSC.W","FITC.A","PE.A","PerCP.Cy5.5.A","PE.Cy7.A","APC.A","APC.H7.A","V450.A","V500.A"]
        df_result = df.ix[df.AccessionNumber==sample,]

        for indx,row in df_result.iterrows():
            df_sample[row["Marker"]] = 0
            for findex in row["flow_index"].split(";"):
                df_sample.ix[int(findex.replace(" ",""))-1,row["Marker"]] = 1


        df_sample.to_csv(result_dir+sample+"_matrix.csv",index=False)


        ### combination analysis
        df_sample["combine"] = df_sample[['PE.ASSC.APOS_NEG', 'PerCP.Cy5.5.AAPC.H7.APOS_NEG','PerCP.Cy5.5.AAPC.H7.ANEG_POS', 'APC.ASSC.APOS_POS',
        'APC.ASSC.APOS_NEG', 'PE.AAPC.APOS_POS', 'PE.AAPC.ANEG_POS']].apply(concatRows,axis=1)
        combination = df_sample.groupby("combine").apply(lambda x : x.shape[0])
        df_combination = pd.DataFrame(combination).reset_index()
        df_combination.columns = ["Combination","Count"]
        df_combination.sort_values("Count",ascending=False,inplace=True)
        df_combination.to_csv(result_dir+sample+"_combination.txt",index=False,sep="\t")


        ### for plotting combination with cd45 only

        selected = ['1010000','1100000','1100110','0000101','1000000','0000001','0001001','1000110','0010000','0100000']


        combination_name = {
        "1010000"    : "CD3+/CD4+/low SSC ~ 3.7%",
        "1100000"    : "CD3+/CD8+/low SSC ~ 2.3%",
        "1100110"    : "CD3+/CD8+/CD56+/low SSC ~ 1.2%",
        "0000101"    : "CD56+/low SSC ~ 1.0%",
        "1000000"    : "CD3+/CD4-/CD8-/low SSC ~ 0.7%",
        "0000001"    : "CD3-/CD56+ ~ 0.6%",
        "0001001"    : "CD3-/CD56+/high SSC ~ 0.6%",
        "1000110"    : "CD3+/CD56+/low SSC ~ 0.4%",
        "0010000"    : "CD3-/CD4+ ~ 0.4%",
        "0100000"    : "CD3-/CD8+ ~ 0.3%"  }

        df_sample2 = df_sample[df_sample["combine"].isin(selected)]

        df_sample2["markers"] = [combination_name[x] for x in df_sample2["combine"].values ]


        g = sns.scatterplot(x="V500.A", y="SSC.A",hue="markers",data=df_sample2,s=10)
        ylabels = ['{:,.0f}'.format(x) + 'K' for x in g.get_yticks()/1000]
        g.set_yticklabels(ylabels)
        plt.xlabel("CD45",fontsize=20)
        plt.ylabel("SSC.A",fontsize=20)
        plt.legend(prop={'size': 16},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+sample+"Tcell_combination_p1.png",bbox_inches = "tight",dpi=300);plt.close()

        g = sns.scatterplot(x="PE.A", y="SSC.A",hue="markers",data=df_sample2,s=10)
        ylabels = ['{:,.0f}'.format(x) + 'K' for x in g.get_yticks()/1000]
        g.set_yticklabels(ylabels)
        plt.xlabel("CD3",fontsize=20)
        plt.ylabel("SSC.A",fontsize=20)
        plt.legend(prop={'size': 16},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+sample+"Tcell_combination_p2.png",bbox_inches = "tight",dpi=300);plt.close()

        sns.scatterplot(x="PerCP.Cy5.5.A", y="APC.H7.A",hue="markers",data=df_sample2,s=10)
        plt.xlabel("CD8",fontsize=20)
        plt.ylabel("CD4",fontsize=20)
        plt.legend(prop={'size': 16},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+sample+"Tcell_combination_p3.png",bbox_inches = "tight",dpi=300);plt.close()

        sns.scatterplot(x="PE.A", y="APC.A",hue="markers",data=df_sample2,s=10)
        plt.xlabel("CD3",fontsize=20)
        plt.ylabel("CD56",fontsize=20)
        plt.legend(prop={'size': 16},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+sample+"Tcell_combination_p4.png",bbox_inches = "tight",dpi=300);plt.close()

def plotBcells(result_file,result_dir):

    df = pd.read_csv(result_file,header=None)
    df.columns =['AccessionNumber','Marker','cell_count','cell_proportion_auto','flow_index']
    df["Marker"] = [x.replace("[","").replace("]","").replace(" ","") for x in df.Marker.values]
    df['AccessionNumber'] = [ x.replace(" ","") for x in df.AccessionNumber.values ]


    cell_marker = {
   'V450.ASSC.APOS_NEGB': 'CD3+SSC-',
   'PE.AFITC.APOS_NEGB-L':  'CD4+CD8-',
   'PE.AFITC.ANEG_POSB-K': 'CD8+CD4-' ,
   'V450.APerCP.Cy5.5.APOS_POSB-CD5' : 'CD56+SSC+',
   'APC.ASSC.APOS_NEGCD10' : 'CD56+SSC-',
   'PE.AFITC.APOS_NEGCD10-L' :  'CD3+SSC-',
   'PE.AFITC.ANEG_POSCD10-K' : 'CD3+SSC-',
   'APC.H7.ASSC.APOS_POSPLASMA': 'CD5'
   }

    for sample in df.AccessionNumber.unique():

        print("processing..."+sample)

        df_sample = pd.read_csv( result_dir+" "+sample+" _trans_mod.csv",header=None)
        df_sample.columns = ["Time","FSC.A","FSC.H","FSC.W","SSC.A","SSC.H","SSC.W","FITC.A","PE.A","PerCP.Cy5.5.A","PE.Cy7.A","APC.A","APC.H7.A","V450.A","V500.A"]
        df_result = df.ix[df.AccessionNumber==sample,]

        for indx,row in df_result.iterrows():
            df_sample[row["Marker"]] = 0
            for findex in row["flow_index"].split(";"):
                df_sample.ix[int(findex.replace(" ",""))-1,row["Marker"]] = 1


        # df_sample.to_csv(result_dir+sample+"_matrix.csv",index=False)

        ### combination analysis
        df_sample["combine"] = df_sample[['V450.ASSC.APOS_NEGB', 'PE.AFITC.APOS_NEGB-L',
           'PE.AFITC.ANEG_POSB-K', 'V450.APerCP.Cy5.5.APOS_POSB-CD5',
           'APC.ASSC.APOS_NEGCD10', 'PE.AFITC.APOS_NEGCD10-L',
           'PE.AFITC.ANEG_POSCD10-K', 'APC.H7.ASSC.APOS_POSPLASMA']].apply(concatRowsBcells,axis=1)
        combination = df_sample.groupby("combine").apply(lambda x : x.shape[0])
        df_combination = pd.DataFrame(combination).reset_index()
        df_combination.columns = ["Combination","Count"]
        df_combination.sort_values("Count",ascending=False,inplace=True)
        df_combination.to_csv(result_dir+sample+"_combination.txt",index=False,sep="\t")

        selected = ['10100000','00001000','00000001','11010000','10001000','00010000','00001001','10001001']

        combination_name = {
        "10100000"    : "CD19+/Kappa/low SSC ~ 0.7%",
        "00001000"    : "CD19-/CD10+/low SSC ~ 0.5%",
        "00000001"    : "CD19-/CD38+ ~ 0.4%",
        "11010000"    : "CD19+/CD5+/Lambda ~ 0.4%",
        "10001000"    : "CD19+/CD10+/low SSC ~ 0.2%",
        "00010000"    : "CD19+/CD5+ ~ 0.2%",
        "00001001"    : "CD19-/CD10+/CD38+/low SSC ~ 0.2%",
        "10001001"    : "CD19+/CD10+/CD38+/low SSC ~ 0.1%"
        }

        df_sample2 = df_sample[df_sample["combine"].isin(selected)]

        df_sample2["markers"] = [combination_name[x] for x in df_sample2["combine"].values ]


        g = sns.scatterplot(x="V500.A", y="SSC.A",hue="markers",data=df_sample2,s=10)
        ylabels = ['{:,.0f}'.format(x) + 'K' for x in g.get_yticks()/1000]
        g.set_yticklabels(ylabels)
        plt.ylabel("SSC.A",fontsize=20)
        plt.xlabel("CD45",fontsize=20)
        plt.legend(prop={'size': 16},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+sample+"Bcell_combination_p1.png",bbox_inches = "tight",dpi=300);plt.close()

        g = sns.scatterplot(x="V450.A", y="SSC.A",hue="markers",data=df_sample2,s=10)
        ylabels = ['{:,.0f}'.format(x) + 'K' for x in g.get_yticks()/1000]
        g.set_yticklabels(ylabels)
        plt.ylabel("SSC.A",fontsize=20)
        plt.xlabel("CD19",fontsize=20)
        plt.legend(prop={'size': 16},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+sample+"Bcell_combination_p2.png",bbox_inches = "tight",dpi=300);plt.close()

        sns.scatterplot(x="PE.A", y="FITC.A",hue="markers",data=df_sample2,s=10)
        plt.xlabel("lambda",fontsize=20)
        plt.ylabel("kappa",fontsize=20)
        plt.legend(prop={'size': 16},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+sample+"Bcell_combination_p3.png",bbox_inches = "tight",dpi=300);plt.close()

def plotCD34cells(result_file,result_dir):

    df = pd.read_csv(result_file,header=None)
    df.columns =['AccessionNumber','Marker','cell_count','cell_proportion_auto','flow_index']
    df["Marker"] = [x.replace("[","").replace("]","").replace(" ","") for x in df.Marker.values]
    df['AccessionNumber'] = [ x.replace(" ","") for x in df.AccessionNumber.values ]


    cell_marker = {
           'APC.ASSC.APOS_NEGCD34': 'CD34+'
           }

    for sample in df.AccessionNumber.unique():

        print("processing..."+sample)

        df_sample = pd.read_csv( result_dir+" "+sample+" _trans_mod.csv",header=None)
        df_sample.columns = ["Time","FSC.A","FSC.H","FSC.W","SSC.A","SSC.H","SSC.W","FITC.A","PE.A","PerCP.Cy5.5.A","PE.Cy7.A","APC.A","APC.H7.A","V450.A","V500.A"]
        df_result = df.ix[df.AccessionNumber==sample,]

        for indx,row in df_result.iterrows():
            df_sample[row["Marker"]] = 0
            for findex in row["flow_index"].split(";"):
                df_sample.ix[int(findex.replace(" ",""))-1,row["Marker"]] = 1

        df_sample.columns = ['Time', 'FSC.A', 'FSC.H', 'FSC.W', 'SSC.A', 'SSC.H', 'SSC.W', 'FITC.A',
        'PE.A', 'PerCP.Cy5.5.A', 'PE.Cy7.A', 'APC.A', 'APC.H7.A', 'V450.A','V500.A', 'markers']

        df_sample["markers"] = [ "CD34+ ~ 1.2%" if x==1 else "CD34-" for x in df_sample.markers.values]

        g = sns.scatterplot(x="V500.A", y="SSC.A",hue="markers",data=df_sample,s=10, palette=['black','red'])
        ylabels = ['{:,.0f}'.format(x) + 'K' for x in g.get_yticks()/1000]
        g.set_yticklabels(ylabels)
        plt.ylabel("SSC.A",fontsize=20)
        plt.xlabel("CD45",fontsize=20)
        plt.legend(prop={'size': 16},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+sample+"CD34cell_1.png",bbox_inches = "tight",dpi=300);plt.close()
        # df_sample.to_csv(result_dir+sample+"_matrix.csv",index=False)


cell_type = sys.argv[1]
result_file = sys.argv[2]
result_dir = sys.argv[3]

print("processing comprehensive expression analysis ")
print("cell type is "+cell_type)
print("result file is "+result_file)

if cell_type == "Tcells":
    plotTcells(result_file,result_dir)
elif cell_type =="Bcells":
    plotBcells(result_file,result_dir)
elif cell_type =="CD34cells":
    plotCD34cells(result_file,result_dir)
