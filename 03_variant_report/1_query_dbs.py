import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import pandas as pd
import seaborn as sns
import numpy as np
import mysql.connector
from mysql.connector import Error
from mysql.connector import errorcode

def getDBInfo(table,chr,pos,ref,alt):
    present = 0
    try:
        connection = mysql.connector.connect(host='localhost',database='', user='', password='')
        cursor = connection.cursor(prepared=True)
        if table=='cosmic':
            query = "select cosmicID from db_cosmic_grch37v88 where chr= %s and pos=%s and ref=%s and alt=%s"
        elif table=='clinvar':
            query = "select clinvarID from db_clinvar_42019 where chr= %s and pos=%s and ref=%s and alt=%s"
        elif table=='gnomad':
            query = "select AF from db_gnomad_r211 where chr= %s and pos=%s and ref=%s and alt=%s"
        elif table=='g1000':
            query = "select g1000ID from db_g1000_phase3v1 where chr= %s and pos=%s and ref=%s and alt=%s"
        cursor.execute(query, (chr,pos,ref,alt))

        for row in cursor:
            res = [el.decode('utf-8') if type(el) is bytearray else el for el in row]
            if len(res) >= 1: present=1

    except mysql.connector.Error as error :
        print("Failed to update record to database:{}".format(error))
        connection.rollback()

    finally:
        if(connection.is_connected()):
            connection.close()

    return present


sample="COLO-829_S5_COLO-829BL_S4.25REFONLY.txt"
# df = pd.read_csv(sample+"_snp_indel_filter.vcf",sep='\t')
df = pd.read_csv(sample,sep='\t')
df_mut = df[['CHROM','POS','REF','ALT']]
# df_mut = df[['chrom','position','ref','var']]
df_mut.columns = ['chrom','position','ref','var']

cosmic_ids=[]
clinvar_ids=[]
gnomad_ids=[]
g1000_ids=[]

for indx, row in df_mut.iterrows():
    print(row['chrom'],row['position'],row['ref'],row['var'])
    cosmic_ids.append(getDBInfo('cosmic',row['chrom'],row['position'],row['ref'],row['var']))
    clinvar_ids.append(getDBInfo('clinvar',row['chrom'],row['position'],row['ref'],row['var']))
    gnomad_ids.append(getDBInfo('gnomad',row['chrom'],row['position'],row['ref'],row['var']))
    g1000_ids.append(getDBInfo('g1000',row['chrom'],row['position'],row['ref'],row['var']))

df_mut['cosmic']=cosmic_ids
df_mut['clinvar']=clinvar_ids
df_mut['gnomad']=gnomad_ids
df_mut['g1000']=gnomad_ids


df_mut.to_csv(sample+".dbout.tsv",sep='\t',index=False)


total = df_mut.shape[0]
all_3 = df_mut[(df_mut['cosmic']==1) & (df_mut['clinvar']==1) & (df_mut['gnomad']==1)].shape[0]
all_none = df_mut[(df_mut['cosmic']==0) & (df_mut['clinvar']==0) & (df_mut['gnomad']==0)].shape[0]
all_2_cosmic_clinvar = df_mut[(df_mut['cosmic']==1) & (df_mut['clinvar']==1) & (df_mut['gnomad']==0)].shape[0]
all_2_cosmic_gnomad = df_mut[(df_mut['cosmic']==1) & (df_mut['clinvar']==0) & (df_mut['gnomad']==1)].shape[0]
all_2_clinvar_gnomad =df_mut[(df_mut['cosmic']==0) & (df_mut['clinvar']==1) & (df_mut['gnomad']==1)].shape[0]
unq_cosmic = df_mut[(df_mut['cosmic']==1) & (df_mut['clinvar']==0) & (df_mut['gnomad']==0)].shape[0]
unq_clinvar = df_mut[(df_mut['cosmic']==0) & (df_mut['clinvar']==1) & (df_mut['gnomad']==0)].shape[0]
unq_gnomad = df_mut[(df_mut['cosmic']==0) & (df_mut['clinvar']==0) & (df_mut['gnomad']==1)].shape[0]
unknown = total - all_3- all_2_cosmic_clinvar - all_2_cosmic_gnomad - all_2_clinvar_gnomad - unq_cosmic - unq_clinvar - unq_gnomad


venn3(subsets = (unq_cosmic,unq_clinvar, all_2_cosmic_clinvar, unq_gnomad,all_2_cosmic_gnomad,all_2_clinvar_gnomad,all_3), set_labels = ('cosmic', 'clinvar', 'gnomad'))
plt.savefig('dbs_venn_bcm')
plt.close()
