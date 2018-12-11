import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib_venn import venn3

import mysql.connector
from mysql.connector import Error
from mysql.connector import errorcode


df = pd.read_csv("COV_7_8_snp_indel_filter.tsv",sep='\t')
df_mut = df[['chrom','position','ref','var']]


def getDBInfo(table,chr,pos,ref,alt):
    present = 0
    try:
        connection = mysql.connector.connect(host='localhost',database='ngs_test', user='', password='')
        cursor = connection.cursor(prepared=True)
        if table=='cosmic':
            query = "select cosmicID from db_cosmic_grch37v86 where chr= %s and pos=%s and ref=%s and alt=%s"
        elif table=='clinvar':
            query = "select clinvarID from db_clinvar where chr= %s and pos=%s and ref=%s and alt=%s"
        elif table=='gnomad':
            query = "select AF from db_gnomad where chr= %s and pos=%s and ref=%s and alt=%s"

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


cosmic_ids=[]
clinvar_ids=[]
gnomad_ids=[]
for indx, row in df_mut.iterrows():
    cosmic_ids.append(getDBInfo('cosmic',row['chrom'],row['position'],row['ref'],row['var']))
    clinvar_ids.append(getDBInfo('clinvar',row['chrom'],row['position'],row['ref'],row['var']))
    gnomad_ids.append(getDBInfo('gnomad',row['chrom'],row['position'],row['ref'],row['var']))

df_mut['cosmic']=cosmic_ids
df_mut['clinvar']=clinvar_ids
df_mut['gnomad']=gnomad_ids

df_mut.to_csv("snp_dbs.tsv",sep='\t',index=False)


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



# Make the diagram
# venn3(subsets = (10, 8, 22, 6,9,4,2))

venn3(subsets = (unq_cosmic,unq_clinvar, all_2_cosmic_clinvar, unq_gnomad,all_2_cosmic_gnomad,all_2_clinvar_gnomad,all_3), set_labels = ('cosmic', 'clinvar', 'gnomad'))
plt.savefig('dbs_venn')
plt.close()
