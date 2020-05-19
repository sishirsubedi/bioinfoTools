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


df = pd.read_excel("data/Final_Combined_Results_v2.xlsx")


#### plot violin


df_violin = pd.melt(df, id_vars =["AccessionNumber","Marker","Pass/Fail","Order"],value_vars =["cell_proportion_manual","cell_proportion_auto"])

df_violin.columns = ['AccessionNumber', 'Marker', 'Pass/Fail', 'Order','method', '% positive cells']
method ={
"cell_proportion_manual":"manual",
"cell_proportion_auto":"auto"
}
df_violin["method"] = [method[x] for x in df_violin["method"].values]
df_violin.sort_values("Order",inplace=True)


plt.figure(figsize=(30,15))
splt = sns.violinplot(x="Marker", y="% positive cells", hue="method",data=df_violin, split=True, palette="Set2",scale="count")
degrees = 25
plt.xticks(rotation=degrees)
splt.set_xticklabels(splt.get_xticklabels(),fontsize=20)
splt.set_xlabel('',fontsize=25)
splt.set_ylabel('% positive cells',fontsize=25)
splt.legend(loc=0, prop={'size': 25})
plt.savefig("Figure1b.png",dpi=300);plt.close()


############### cv bar plot

df_manual = df_manual[df_manual["Pass/Fail"]=="P"]


df_manual_melt = pd.melt(df, id_vars =["AccessionNumber","Marker","Pass/Fail"],value_vars =["cell_proportion_manual","cell_proportion_auto"])

method ={
"cell_proportion_manual":"manual",
"cell_proportion_auto":"auto"
}
df_manual_melt["Method"] = [method[x] for x in df_manual_melt["variable"].values]

df_cv = df_manual_melt.groupby(["Marker","Method"]).agg({'value': ['mean', 'std',]})



########### pass fail stack bar

df.sort_values("Order",inplace=True)
t = df.groupby(["Marker","Order","Pass/Fail"]).count()["AccessionNumber"].unstack()
t.reset_index(inplace=True)
t.sort_values("Order",ascending=False,inplace=True)
t = t[["Marker",'P','F']]
t.columns = ["Marker",'Pass','Fail']
t.set_index("Marker",inplace=True)

colors = ['seagreen', 'royalblue']
t.plot(kind='barh', stacked=True,color=colors,figsize=(30,20))
ax = plt.gca()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Total cases',fontsize=25)
plt.ylabel('')
plt.legend(prop={'size': 25},loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
plt.savefig("Figure2a.png",dpi=300);plt.close()




####### corrrelation table

df.sort_values("Order",inplace=True)
df = df[df["Pass/Fail"]=="P"]
for marker in df.Marker.unique():

    df_sample = df[df["Marker"]==marker]
    r_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[0]
    p_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[1]

    print(marker + " , pearson-r : "+str(round(r_val,4)))


##### correlation plot best 6

df.sort_values("Order",inplace=True)
df = df[df["Pass/Fail"]=="P"]


subplots =[
'CD3+/low SSC',
'CD4+/CD8-',
'CD8+/CD4-',
'CD3-/CD56+ (low SSC)',
'CD3+/CD56+ (low SSC)',
'CD19+/low SSC',
'Kappa (CD19+)',
'Lambda (CD19+)',
'CD10+/low SSC',
'CD34+/low SSC'
]

# plt.rcParams.update({'font.size': 12})
# fig, axes = plt.subplots(3,3, sharex=False, sharey=False,figsize=(16,10))
# for i, ax in zip(range(9), axes.flat):
#         marker = subplots[i]
#
#         # df_sample = df[( (df["Marker"]==marker) & ( df['Category']=='p' ))]
#
#         df_sample = df[ (df["Marker"]==marker) ]
#
#
#         r_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[0]
#         p_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[1]
#
#         sns.scatterplot(x="cell_proportion_manual", y="cell_proportion_auto",
#         hue="Marker",data=df_sample,legend=False,palette=['black'],ax=ax)
#         plot_text = "Spearman r = "+str(round(r_val,4))+ "\n n="+str(df_sample.shape[0])+" pairs"
#         ax.set_ylabel(" % of "+marker+" cells (flowDensity)")
#         ax.set_xlabel(" % of "+marker+" cells (manual)")
#         ax.text(0.3, 0.9,plot_text, ha='center', va='center', transform=ax.transAxes)
#
#
# fig.tight_layout(pad=2.0)
# plt.savefig("Figure2b.png");plt.close()


plt.rcParams.update({'font.size': 12})
fig, axes = plt.subplots(2,5, sharex=False, sharey=False,figsize=(16,8))
for i, ax in zip(range(10), axes.flat):
        marker = subplots[i]

        # df_sample = df[( (df["Marker"]==marker) & ( df['Category']=='p' ))]

        df_sample = df[ (df["Marker"]==marker) ]


        r_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[0]
        p_val = stats.spearmanr(df_sample["cell_proportion_manual"].values,df_sample["cell_proportion_auto"].values)[1]

        sns.scatterplot(x="cell_proportion_manual", y="cell_proportion_auto",
        hue="Marker",data=df_sample,legend=False,palette=['black'],ax=ax)
        plot_text = "r:"+str(round(r_val,4))+ "\n n="+str(df_sample.shape[0])+" pairs"
        # ax.set_ylabel(" % of "+marker+" cells (flowDensity)")
        ax.set_xlabel(marker)
        ax.set_ylabel("")
        ax.text(0.3, 0.9,plot_text, ha='center', va='center', transform=ax.transAxes)


fig.text(0.05, 0.0, 'manual->', ha='center', va='center')
fig.text(0.0, 0.08, 'flowDensity->', ha='center', va='center', rotation='vertical')
fig.tight_layout(pad=0.5)
plt.savefig("Figure2b.png",dpi=300, bbox_inches = "tight");plt.close()
