
import pandas as pd

file="output.snvs.vcf.filter.txt"
# file ="output.indels.vcf.filter.txt"
type=file.split('.')[1]



def checkAlleleFreq(counts,th):
    arr = [int(x.split(',')[0])+int(x.split(',')[1]) for x in counts]
    max_indx = arr.index(max(arr))
    check = 0
    for i in range(len(arr)):
        if (i != max_indx) and ( arr[i] > th):
            check = 1
            break
    return check


def parseDF(df,type):
    # df=df[((df['FILTER'] == "PASS") | (df['FILTER'] =="LowEVS"))]
    # df=df[(df['FILTER'] == "PASS")]
    df=df[((df['NORMAL.DP']>=20) & (df['TUMOR.DP']>=20) & (df['SomaticEVS']>=4 ) ) ]
    if type =='snvs':
        filter_indx=[]
        for i,r in df.iterrows():
            if (checkAlleleFreq([r['NORMAL.AU'],r['NORMAL.CU'],r['NORMAL.GU'],r['NORMAL.TU']],0)):
                filter_indx.append(i)
        df.drop(filter_indx,inplace=True)
    return df

df = pd.read_csv(file,sep='\t')
df_res = parseDF(df,type)
df_res.to_csv(type+".strelka.filter.txt",sep='\t',index=False)
