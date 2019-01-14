import sys
import pandas as pd

def checkPosition(x1,y1,x2,y2,olap):
    x1 += olap
    y1 -= olap
    if (x2 >= x1 and x2 <= y1) or (y2 >= x1 and y2 <= y1) :
        return True
    else:
        return False

def filterExons(df_e,df_d,olap):
    filter=[]
    for indx,row in df_d.iterrows():
        df_e_mini = df_e[(df_e['chr']==row['chr']) &
        ( ( row['start'] >= (df_e['start']+olap)) & ( row['start'] <= (df_e['end']-olap) ) ) |
        ( ( row['end'] >= (df_e['start']+olap)) & ( row['end'] <= (df_e['end']-olap) ) ) ]
        if df_e_mini.shape[0]>=1:
            print([row['chr'],row['start'],row['end'],df_e_mini.shape[0]])
            filter.append([row['chr'],row['start'],row['end'],df_e_mini.shape[0]])
    return filter

exons= "/home/hhadmin/exome_pipeline/01_bamQC/ucsc_hg19_exons_only_coordinates.txt"
design= "/home/hhadmin/exome_pipeline/01_bamQC/cre_design.bed"
output="/home/hhadmin/exome_pipeline/01_bamQC/"

df_exons= pd.read_csv(exons,sep='\t',header=None)
df_exons=df_exons.iloc[:,0:3]
df_exons.columns=["chr","start",'end']
df_exons.set_index(['chr','start','end'])

df_design= pd.read_csv(design,sep='\t',header=None)
df_design.columns=["chr","start",'end']
df_design.set_index(['chr','start','end'])


result= filterExons(df_exons,df_design,2)
pd.DataFrame(result).to_csv(output+"cre_design_exonOnly_withexonCount.bed",sep='\t',index=False)
