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


def plotTcells():
    result_file = "/net/fs01.cluster.com/home/hhadmin/flowCyto/indepen_analysis/results/t_cells_result/ T_CELLS_result.csv"
    result_dir = "/net/fs01.cluster.com/home/hhadmin/flowCyto/indepen_analysis/results/t_cells_result/"

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

        #if sample == "FLW1919":break

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
        "1010000"    : "CD3+/CD4+/low SSC",
        "1100000"    : "CD3+/CD8+/low SSC",
        "1100110"    : "CD3+/CD8+/CD56+/low SSC",
        "0000101"    : "CD56+/low SSC",
        "1000000"    : "CD3+/CD4-/CD8-/low SSC",
        "0000001"    : "CD3-/CD56+",
        "0001001"    : "CD3-/CD56+/high SSC",
        "1000110"    : "CD3+/CD56+/low SSC",
        "0010000"    : "CD3-/CD4+",
        "0100000"    : "CD3-/CD8+"  }

        df_sample2 = df_sample[df_sample["combine"].isin(selected)]

        df_sample2["markers"] = [combination_name[x] for x in df_sample2["combine"].values ]


        sns.scatterplot(x="V500.A", y="SSC.A",hue="markers",data=df_sample2,s=10)
        plt.xlabel("CD45")
        plt.legend(prop={'size': 12},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+"Tcell_combination_p1.png",bbox_inches = "tight",dpi=300);plt.close()

        sns.scatterplot(x="PE.A", y="SSC.A",hue="markers",data=df_sample2,s=10)
        plt.xlabel("CD3")
        plt.legend(prop={'size': 12},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+"Tcell_combination_p2.png",bbox_inches = "tight",dpi=300);plt.close()

        sns.scatterplot(x="PerCP.Cy5.5.A", y="APC.H7.A",hue="markers",data=df_sample2,s=10)
        plt.xlabel("CD8")
        plt.ylabel("CD4")
        plt.legend(prop={'size': 12},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+"Tcell_combination_p3.png",bbox_inches = "tight",dpi=300);plt.close()

        sns.scatterplot(x="PE.A", y="APC.A",hue="markers",data=df_sample2,s=10)
        plt.xlabel("CD3")
        plt.ylabel("CD56")
        plt.legend(prop={'size': 12},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+"Tcell_combination_p4.png",bbox_inches = "tight",dpi=300);plt.close()



        #### for plotting

        df_sample["CD3_Marker"]=""
        df_sample.ix[df_sample["PE.ASSC.APOS_NEG"]==1,["CD3_Marker"]]="CD3"
        df_sample.ix[df_sample["PerCP.Cy5.5.AAPC.H7.APOS_NEG"]==1,["CD3_Marker"]]="CD8"
        df_sample.ix[df_sample["PerCP.Cy5.5.AAPC.H7.ANEG_POS"]==1,["CD3_Marker"]]="CD4"


        df_sample["CD56_Marker"]=""
        df_sample.ix[df_sample["APC.ASSC.APOS_NEG"]==1,["CD56_Marker"]]="CD56"
        df_sample.ix[df_sample["PE.AAPC.APOS_POS"]==1,["CD56_Marker"]]="TLGL"
        df_sample.ix[df_sample["PE.AAPC.ANEG_POS"]==1,["CD56_Marker"]]="NK"


        sns.scatterplot(x="PE.A", y="SSC.A",hue="CD3_Marker",data=df_sample,s=10)
        plt.savefig("Tcell.png");plt.close()

        sns.scatterplot(x="PerCP.Cy5.5.A", y="APC.H7.A",hue="CD3_Marker",data=df_sample,s=10)
        plt.savefig("Tcell_2.png");plt.close()


        sns.scatterplot(x="V500.A", y="SSC.A",hue="CD3_Marker",data=df_sample[df_sample["PE.ASSC.APOS_NEG"]==1],s=10)
        plt.savefig("Tcell_3.png");plt.close()

        sns.scatterplot(x="V500.A", y="SSC.A",hue="CD3_Marker",data=df_sample,s=10)
        plt.savefig("Tcell_4.png");plt.close()


        sns.scatterplot(x="PE.A", y="APC.A",hue="CD56_Marker",data=df_sample[df_sample["APC.ASSC.APOS_NEG"]==1],s=10)
        plt.savefig("Tcell_5.png");plt.close()


        sns.scatterplot(x="PE.A", y="APC.A",hue="CD56_Marker",data=df_sample,s=10)
        plt.savefig("Tcell_6.png");plt.close()


        ### cd3 to all types

        df_sample=df_sample[df_sample["PE.ASSC.APOS_NEG"]==1]

def plotBcells():

    result_file = "/home/hhadmin/flowCyto/indepen_analysis/results/b_cells_result/ B_CELLS_result.csv"
    result_dir = "/home/hhadmin/flowCyto/indepen_analysis/results/b_cells_result/"

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

        #if sample == "FLW1919":break

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
        "10100000"    : "CD19+/Kappa/low SSC",
        "00001000"    : "CD19-/CD10+/low SSC",
        "00000001"    : "CD19-/CD38+",
        "11010000"    : "CD19+/CD5+/Lambda",
        "10001000"    : "CD19+/CD10+/low SSC",
        "00010000"    : "CD19+/CD5+",
        "00001001"    : "CD19-/CD10+/CD38+/low SSC",
        "10001001"    : "CD19+/CD10+/CD38+/low SSC"
        }

        df_sample2 = df_sample[df_sample["combine"].isin(selected)]

        df_sample2["markers"] = [combination_name[x] for x in df_sample2["combine"].values ]


        sns.scatterplot(x="V500.A", y="SSC.A",hue="markers",data=df_sample2,s=10)
        plt.xlabel("CD45")
        plt.legend(prop={'size': 12},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+"Bcell_combination_p1.png",bbox_inches = "tight",dpi=300);plt.close()

        sns.scatterplot(x="V450.A", y="SSC.A",hue="markers",data=df_sample2,s=10)
        plt.xlabel("CD19")
        plt.legend(prop={'size': 12},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+"Bcell_combination_p2.png",bbox_inches = "tight",dpi=300);plt.close()

        sns.scatterplot(x="PE.A", y="FITC.A",hue="markers",data=df_sample2,s=10)
        plt.xlabel("lambda")
        plt.ylabel("kappa")
        plt.legend(prop={'size': 12},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
        plt.savefig(result_dir+"Bcell_combination_p3.png",bbox_inches = "tight",dpi=300);plt.close()


def plotcd34():

    result_file = "/home/hhadmin/flowCyto/indepen_analysis/results/cd34_cells_result/ CD34_CELLS_result.csv"
    result_dir = "/home/hhadmin/flowCyto/indepen_analysis/results/cd34_cells_result/"

    df = pd.read_csv(result_file,header=None)
    df.columns =['AccessionNumber','Marker','cell_count','cell_proportion_auto','flow_index']
    df["Marker"] = [x.replace("[","").replace("]","").replace(" ","") for x in df.Marker.values]
    df['AccessionNumber'] = [ x.replace(" ","") for x in df.AccessionNumber.values ]


    cell_marker = {
           'APC.ASSC.APOS_NEGCD34': 'CD34+'
           }


    sample = "FLW195056"

    print("processing..."+sample)

    df_sample = pd.read_csv( result_dir+" "+sample+" _trans_mod.csv",header=None)
    df_sample.columns = ["Time","FSC.A","FSC.H","FSC.W","SSC.A","SSC.H","SSC.W","FITC.A","PE.A","PerCP.Cy5.5.A","PE.Cy7.A","APC.A","APC.H7.A","V450.A","V500.A"]
    df_result = df.ix[df.AccessionNumber==sample,]

    for indx,row in df_result.iterrows():
        df_sample[row["Marker"]] = 0
        for findex in row["flow_index"].split(";"):
            df_sample.ix[int(findex.replace(" ",""))-1,row["Marker"]] = 1

    df_sample.columns = ['Time', 'FSC.A', 'FSC.H', 'FSC.W', 'SSC.A', 'SSC.H', 'SSC.W', 'FITC.A',
       'PE.A', 'PerCP.Cy5.5.A', 'PE.Cy7.A', 'APC.A', 'APC.H7.A', 'V450.A',
       'V500.A', 'markers']

    df_sample["markers"] = [ "CD34+" if x==1 else "CD34-" for x in df_sample.markers.values]

    sns.scatterplot(x="V500.A", y="SSC.A",hue="markers",data=df_sample,s=10)
    plt.xlabel("CD45")
    plt.legend(prop={'size': 12},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=False)
    plt.savefig(result_dir+"CD34cell_1.png",bbox_inches = "tight",dpi=300);plt.close()
    # df_sample.to_csv(result_dir+sample+"_matrix.csv",index=False)
