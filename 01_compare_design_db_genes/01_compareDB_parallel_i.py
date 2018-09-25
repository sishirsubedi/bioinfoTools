
import pandas as pd
import optparse
import sys
import threading
from multiprocessing import Process, Queue


def compareDB(database,db_rows,bedfile, bedfile_rows):

    df_db = pd.read_csv(database, sep='\t', header=0, skiprows=db_rows)
    # print(len(df_db.index))
    filter_indx = df_db[df_db['CHROM'].str.contains("MT|NW_00")].index
    # print (len(filter_indx))
    df_db.drop(filter_indx, inplace=True)
    # print(len(df_db.index))
    # print(df_db.tail())

    df_bedfile = pd.read_csv(bedfile, sep='\t', header=None, skiprows=bedfile_rows)
    df_bedfile = df_bedfile.iloc[:, 0:3]
    df_bedfile.columns =['file_CHROM','file_START','file_STOP']

    def find_match(chromosome, df_chromosome, out_queue):

        for indx,row in df_chromosome.iterrows():
            coverage=0
            temp=[]
            if len(df_bedfile[(df_bedfile['file_START'] <= row['POS']) & (df_bedfile['file_STOP'] >= row['POS']) & (df_bedfile['file_CHROM'] == row['CHROM'])].index)==1:
                coverage = 1
            temp = [x for x in row.values]
            temp.append(coverage)
            out_queue.put(temp)



    chromosomes = df_db['CHROM'].unique()
    counter = 0
    combine = []
    totalchrs = len(chromosomes)

    while counter < totalchrs:

        print(counter, totalchrs)

        out_queue = Queue()
        processes = []

        for i in range(3):
            process = Process(target=find_match,
                              args=(chromosomes[counter], df_db[df_db['CHROM'] == chromosomes[counter]], out_queue))
            processes.append(process)
            process.start()
            counter = counter + 1

        for p in processes:
            p.join()

        for i in range(out_queue.qsize()):
            temp = []
            temp = out_queue.get()
            combine.append(temp)



    df_combine = pd.DataFrame(combine)
    df_combine.columns = ['CHROM', 'POS', 'ID ','GENE','FILE_Coverage']
    df_combine.to_csv("out_db_comparison.csv")


try:
    parser = optparse.OptionParser()
    parser.add_option('-d', '--database', help = ' database')
    parser.add_option('-r', '--rows', help = 'rows to skip')
    parser.add_option('-i', '--infile', help = ' input bed file')
    parser.add_option('-n', '--infilerows', help=' input bed file skip rows')
    options,args = parser.parse_args()
    database = options.database
    rows = options.rows
    bedfile = options.infile
    bedfile_rows = options.infilerows
    compareDB(database,rows,bedfile, bedfile_rows)

except TypeError:
	print ("python 02_1_compareDB.py -help for help")