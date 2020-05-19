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

result_file = "data/B_CELLS_RESULT_OUTPUT_BCells_parents/ B_CELLS_result.csv"
result_dir =  "data/B_CELLS_RESULT_OUTPUT_BCells_parents/ "


df = pd.read_csv(result_file,header=None)
df.columns =['AccessionNumber','Marker','cell_count','cell_proportion_auto','flow_index']
df["Marker"] = [x.replace("[","").replace("]","").replace(" ","") for x in df.Marker.values]
df['AccessionNumber'] = [ x.replace(" ","") for x in df.AccessionNumber.values ]

df = df[['AccessionNumber', 'Marker', 'cell_proportion_auto']]
# df.to_csv(result_dir+"B_Cells_parents_compare_result.csv",index=False)
df.to_csv(result_dir+"CD34_compare_result.csv",index=False)

################## manual comparison


df_manual = pd.read_excel("data/Manual_103NormalCases.xlsx")
df_manual.columns = [x.replace(" ","") for x in df_manual.columns]

df_manual_melt = pd.melt(df_manual,id_vars=['AccessionNumber'],value_vars=['PE.ASSC.APOS_NEG',
    'PerCP.Cy5.5.AAPC.H7.ANEG_POS','PerCP.Cy5.5.AAPC.H7.APOS_NEG', 'APC.ASSC.APOS_NEG',
    'PE.AAPC.ANEG_POS','PE.AAPC.APOS_POS', 'APC.ASSC.APOS_POS'])
df_manual_melt.sort_values(by='AccessionNumber',inplace=True)
df_manual_melt.columns = ['AccessionNumber','Marker','cell_proportion_manual']

dfjoin = pd.merge(df_manual_melt,df,on=['AccessionNumber','Marker'],how='outer',indicator='True')
dfjoin.to_csv("T_Cells_compare_result.csv")

sns.scatterplot(x="cell_proportion_manual", y="cell_proportion_auto",hue="Marker",data=dfjoin)
plt.savefig("temp.png");plt.close()

cell_marker={
'PE.ASSC.APOS_NEG' : 'CD3',
'PerCP.Cy5.5.AAPC.H7.APOS_NEG' : 'CD8',
'PerCP.Cy5.5.AAPC.H7.ANEG_POS' : 'CD4',
'APC.ASSC.APOS_POS' : 'CD56 Other',
'APC.ASSC.APOS_NEG' : 'CD56 Lymphs'  ,
'PE.AAPC.APOS_POS' :  'NK Like T LGL',
'PE.AAPC.ANEG_POS' : 'NK'
}

for marker in dfjoin.Marker.unique():

    df_sample = dfjoin[dfjoin["Marker"]==marker]
    r_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[0]
    p_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[1]

    sns.scatterplot(x="cell_proportion_manual", y="cell_proportion_auto",hue="Marker",data=df_sample,legend=False)

    sns.regplot(x="cell_proportion_manual", y="cell_proportion_auto",data=df_sample)



    plt.legend( loc='upper left', labels=["pearson-r : "+str(round(r_val,4)),cell_marker[marker]])
    plt.savefig(marker+"_correlation_plot.png");plt.close()

    print(cell_marker[marker] + " , pearson-r : "+str(round(r_val,4)))



##### expert analysis

cell_marker={
'PE.ASSC.APOS_NEG' : 'CD3',
'PerCP.Cy5.5.AAPC.H7.APOS_NEG' : 'CD8',
'PerCP.Cy5.5.AAPC.H7.ANEG_POS' : 'CD4',
'APC.ASSC.APOS_POS' : 'CD56 Other',
'APC.ASSC.APOS_NEG' : 'CD56 Lymphs'  ,
'PE.AAPC.APOS_POS' :  'NK Like T LGL',
'PE.AAPC.ANEG_POS' : 'NK'
}


df=pd.read_excel("T_Cells_compare_result.xls")
df['cell_type']=""
df.columns = ['AccessionNumber', 'Marker', 'cell_proportion_manual',
       'cell_proportion_auto', 'TRUE', 'Category', 'Note', 'cell_type']

for m in cell_marker.keys():
    df.ix[df["Marker"]==m,'cell_type'] = cell_marker[m]

df.ix[df['cell_type']=="CD56 Lymphs",'cell_type'] = "CD56-Ly"
df.ix[df['cell_type']=="NK Like T LGL",'cell_type'] = "NK-TLGL"
df.ix[df['cell_type']=="CD56 Other",'cell_type'] = "CD56-Ot"


plt.rcParams.update({'font.size': 20})
colors = ['forestgreen', 'firebrick','cornflowerblue']
fig_size = (15,12)
t = df.groupby(['cell_type','Category']).count()['AccessionNumber'].unstack()
t=t.reindex(['CD3','CD4', 'CD8', 'CD56-Ly',  'NK','NK-TLGL','CD56-Ot'])
t=t[['p','f','s']]
t.columns = ['Pass','Failed','Satis']
t.plot(kind='bar', stacked=True,color=colors,figsize=(15,10))
plt.xticks(rotation=45,ha='right')
plt.xlabel("")
plt.ylabel("Number of samples (n=103)")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, fancybox=True, shadow=True)
plt.savefig("Figure1_all_bar_plot.png");plt.close()



plt.rcParams.update({'font.size': 12})

fig, axes = plt.subplots(2,3, sharex=False, sharey=False,figsize=(16,10))

subplots =[
'PE.ASSC.APOS_NEG',
'PerCP.Cy5.5.AAPC.H7.ANEG_POS',
'PerCP.Cy5.5.AAPC.H7.APOS_NEG',
'APC.ASSC.APOS_NEG',
'PE.AAPC.APOS_POS',
'PE.AAPC.ANEG_POS']

for i, ax in zip(range(6), axes.flat):
        marker = subplots[i]

        # df_sample = df[( (df["Marker"]==marker) & ( df['Category']=='p' ))]

        df_sample = df[ (df["Marker"]==marker) ]


        r_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[0]
        p_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[1]

        sns.scatterplot(x="cell_proportion_manual", y="cell_proportion_auto",
        hue="Marker",data=df_sample,legend=False,palette=['black'],ax=ax)
        plot_text = "Spearman r = "+str(round(r_val,4))+ "\n n="+str(df_sample.shape[0])+" pairs"
        ax.set_ylabel(" % of "+cell_marker[marker]+"+ cells (flowDensity)")
        ax.set_xlabel(" % of "+cell_marker[marker]+"+ cells (manual)")
        ax.text(0.3, 0.9,plot_text, ha='center', va='center', transform=ax.transAxes)


fig.tight_layout(pad=2.0)
plt.savefig("Figure2_All_correlation_plot_raw.png");plt.close()




######################## b B_CELLS_result


df_manual = pd.read_excel("data/Manual_103NormalCases_BCells.xlsx")
df_manual.columns = [x.replace(" ","") for x in df_manual.columns]

df_manual_melt = pd.melt(df_manual,id_vars=['AccessionNumber'],value_vars=['V450.ASSC.APOS_NEG',
    'APC.ASSC.APOS_NEG','APC.H7.ASSC.APOS_NEG'])
df_manual_melt.sort_values(by='AccessionNumber',inplace=True)
df_manual_melt.columns = ['AccessionNumber','Marker','cell_proportion_manual']


df= df[df['AccessionNumber'] != 'FLW194995']

dfjoin = pd.merge(df_manual_melt,df,on=['AccessionNumber','Marker'],how='outer',indicator='True')
dfjoin.to_csv("data/B_Cells_parents_compare_result.csv")

sns.scatterplot(x="cell_proportion_manual", y="cell_proportion_auto",hue="Marker",data=dfjoin)
plt.savefig("temp.png");plt.close()

cell_marker={
'V450.ASSC.APOS_NEG' : 'B',
'APC.ASSC.APOS_NEG' : 'CD10',
'APC.H7.ASSC.APOS_NEG' : 'CD38',
}

for marker in dfjoin.Marker.unique():

    df_sample = dfjoin[dfjoin["Marker"]==marker]
    r_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[0]
    p_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[1]

    sns.scatterplot(x="cell_proportion_manual", y="cell_proportion_auto",hue="Marker",data=df_sample,legend=False)

    # sns.regplot(x="cell_proportion_manual", y="cell_proportion_auto",data=df_sample)



    plt.legend( loc='upper left', labels=["pearson-r : "+str(round(r_val,4)),cell_marker[marker]])

    plt.title(cell_marker[marker])

    plt.savefig("data/"+marker+"_correlation_plot.png");plt.close()

    print(cell_marker[marker] + " , pearson-r : "+str(round(r_val,4)))



##### expert analysis

cell_marker={
'PE.ASSC.APOS_NEG' : 'CD3',
'PerCP.Cy5.5.AAPC.H7.APOS_NEG' : 'CD8',
'PerCP.Cy5.5.AAPC.H7.ANEG_POS' : 'CD4',
'APC.ASSC.APOS_POS' : 'CD56 Other',
'APC.ASSC.APOS_NEG' : 'CD56 Lymphs'  ,
'PE.AAPC.APOS_POS' :  'NK Like T LGL',
'PE.AAPC.ANEG_POS' : 'NK'
}


df=pd.read_excel("T_Cells_compare_result.xls")
df['cell_type']=""
df.columns = ['AccessionNumber', 'Marker', 'cell_proportion_manual',
       'cell_proportion_auto', 'TRUE', 'Category', 'Note', 'cell_type']

for m in cell_marker.keys():
    df.ix[df["Marker"]==m,'cell_type'] = cell_marker[m]

df.ix[df['cell_type']=="CD56 Lymphs",'cell_type'] = "CD56-Ly"
df.ix[df['cell_type']=="NK Like T LGL",'cell_type'] = "NK-TLGL"
df.ix[df['cell_type']=="CD56 Other",'cell_type'] = "CD56-Ot"


plt.rcParams.update({'font.size': 20})
colors = ['forestgreen', 'firebrick','cornflowerblue']
fig_size = (15,12)
t = df.groupby(['cell_type','Category']).count()['AccessionNumber'].unstack()
t=t.reindex(['CD3','CD4', 'CD8', 'CD56-Ly',  'NK','NK-TLGL','CD56-Ot'])
t=t[['p','f','s']]
t.columns = ['Pass','Failed','Satis']
t.plot(kind='bar', stacked=True,color=colors,figsize=(15,10))
plt.xticks(rotation=45,ha='right')
plt.xlabel("")
plt.ylabel("Number of samples (n=103)")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, fancybox=True, shadow=True)
plt.savefig("Figure1_all_bar_plot.png");plt.close()



plt.rcParams.update({'font.size': 12})

fig, axes = plt.subplots(2,3, sharex=False, sharey=False,figsize=(16,10))

subplots =[
'PE.ASSC.APOS_NEG',
'PerCP.Cy5.5.AAPC.H7.ANEG_POS',
'PerCP.Cy5.5.AAPC.H7.APOS_NEG',
'APC.ASSC.APOS_NEG',
'PE.AAPC.APOS_POS',
'PE.AAPC.ANEG_POS']

for i, ax in zip(range(6), axes.flat):
        marker = subplots[i]

        # df_sample = df[( (df["Marker"]==marker) & ( df['Category']=='p' ))]

        df_sample = df[ (df["Marker"]==marker) ]


        r_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[0]
        p_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[1]

        sns.scatterplot(x="cell_proportion_manual", y="cell_proportion_auto",
        hue="Marker",data=df_sample,legend=False,palette=['black'],ax=ax)
        plot_text = "Spearman r = "+str(round(r_val,4))+ "\n n="+str(df_sample.shape[0])+" pairs"
        ax.set_ylabel(" % of "+cell_marker[marker]+"+ cells (flowDensity)")
        ax.set_xlabel(" % of "+cell_marker[marker]+"+ cells (manual)")
        ax.text(0.3, 0.9,plot_text, ha='center', va='center', transform=ax.transAxes)


fig.tight_layout(pad=2.0)
plt.savefig("Figure2_All_correlation_plot_raw.png");plt.close()
