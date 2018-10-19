import pandas as pd

outfile=''
infile=''

chunksize = 10 ** 6
with open(outfile,"w") as w:
    header = ['gnomadID','chr','pos','ref','alt','AF']
    w.writelines('\t'.join(str(j) for j in header) + '\n')
    count=0
    for df in pd.read_csv(infile,sep='\t', chunksize=chunksize):
        df[df['filter']=='None']
        df = df[['gnomad-id', 'chr', 'pos', 'ref', 'alt','AF']]
        df = df[df['AF']!='na']
        filter=[]
        for indx,row in df.iterrows():
            alts = row['alt'].split(',')
            if len(alts) == 1:
                alt=alts[0].replace(']','').replace('[','').replace("'",'')
                filter.append([row['gnomad-id'],'chr'+str(row['chr']),row['pos'],row['ref'],alt,row['AF']])
            else:
                af_list=row['AF'].replace(' ','').replace('(','').replace(')','').split(',')
                for i in range(len(alts)):
                    alt=alts[i].replace(']','').replace('[','').replace("'",'')
                    if af_list[i] != '0.0':
                        filter.append([row['gnomad-id'],'chr'+str(row['chr']),row['pos'],row['ref'],alt,af_list[i]])
        w.writelines('\t'.join(str(j) for j in i) + '\n' for i in filter)
        count +=1
        print('finished -',count)
