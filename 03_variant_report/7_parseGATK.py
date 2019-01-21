
import pandas as pd

def calcSum(row):
    vals= row.split(',')
    return int(vals[0])+int(vals[1])

df_main = pd.DataFrame()
cols=[]
# for chr in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']:
for chr in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']:
    print(chr)
    df=pd.read_csv("chr"+chr+"_gatk.output.vcf.filter.txt",sep='\t')
    print(df.head(2))
    df['NORMAL.DP'] = df['NORMAL.AD'].apply(calcSum)
    df['TUMOR.DP'] = df['TUMOR.AD'].apply(calcSum)
    df = df[(df['NORMAL.DP']>=20) & (df['TUMOR.DP']>=20 )]
    df_main = df_main.append(df,ignore_index=True)
    if chr =='1':
        cols=df.columns

print(cols)
df_main.columns=cols
df_main.to_csv("01_gatk_filter_combine.vcf.txt",sep='\t',index=False)
