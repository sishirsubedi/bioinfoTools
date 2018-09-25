
import pandas as pd
import optparse
import sys
import threading
from multiprocessing import Process, Queue




df_clinvar = pd.read_csv(sys.argv[1], sep='\t', skiprows=27)
df_clinvar = pd.read_csv("db_clinvar_mod.vcf", sep='\t', header=0, skiprows=27)

print(len(df_clinvar.index))
filter_indx = df_clinvar[df_clinvar['CHROM'].str.contains("MT|NW_00")].index
print (len(filter_indx))
df_clinvar.drop(filter_indx, inplace=True)
print(len(df_clinvar.index))
print(df_clinvar.tail())


#df_bedfile = pd.read_csv(sys.argv[2], sep='\t', header=None, skiprows=2)
df_bedfile = pd.read_csv("v7_S31285117_Regions.bed", sep='\t', header=None, skiprows=2)
df_bedfile = df_bedfile.iloc[:,0:3]
df_bedfile.columns =['v7_CHROM','v7_START','v7_STOP']


#df_bedfile2 = pd.read_csv(sys.argv[3], sep='\t', header=None, skiprows=2)
df_bedfile2 = pd.read_csv("cre2_S30409818_Regions.bed", sep='\t', header=None, skiprows=2)
df_bedfile2 = df_bedfile2.iloc[:,0:3]
df_bedfile2.columns =['cre2_CHROM','cre2_START','cre2_STOP']


class design:
	def __init__(self, cov, cov_num):
		self.cov=cov
		self.cov_num=cov_num


def find_coverage(df,chromosome,clinvar_position,design,fn):
	match = df[ (df[fn+'_START'] <= clinvar_position) & (df[fn+'_STOP'] >= clinvar_position) & (df[fn+'_CHROM'] == chromosome)]
	if (len(match.index)>=1):
		design.cov = 1
		design.cov_num += 1

def find_match(gene, df_gene,out_queue):

	## no neeed to create a global variable because this variable is read only
	# global df_bedfile
	# global df_bedfile2

	if gene=="nan":return
	if(len(df_gene.index)==0):return

	clinvar_positions = df_gene['POS'].values
	chromosome= df_gene['CHROM'].values[0]

	threads =[]
	d2=design(0,0)
	d1=design(0,0)

	for clinvar_position in clinvar_positions:

		thread1 = threading.Thread(target=find_coverage, args=(df_bedfile,chromosome,clinvar_position,d1,'v7'))
		threads.append(thread1)
		thread1.start()

		thread2 = threading.Thread(target=find_coverage, args=(df_bedfile2,chromosome,clinvar_position,d2,'cre2'))
		threads.append(thread2)
		thread2.start()

	# wait for all threads to finish

	for t in threads:
		t.join()

	out_queue.put([chromosome, gene,len(clinvar_positions), d1.cov, d1.cov_num, d2.cov, d2.cov_num ])



genes=df_clinvar['GENE'].unique()
counter =0
combine=[]

totalgenes=len(genes)

while counter < totalgenes:

	out_queue=Queue()
	processes=[]

	for i in range(2):
		process = Process(target=find_match, args=(genes[counter],df_clinvar[df_clinvar['GENE']==genes[counter]],out_queue))
		processes.append(process)
		process.start()
		counter = counter + 1

	for p in processes:
		p.join()

	for i in range(out_queue.qsize()):
		temp=[]
		temp=out_queue.get()
		combine.append(temp)

	print(counter,totalgenes)

df_combine = pd.DataFrame(combine)

df_combine.columns = ['chromosome', 'gene', 'clinvar_exons','v7','v7_exon', 'cre2', 'cre2_exon']
df_combine.to_csv("out_clinvar_comparison2_cre_v7_bedfiles.csv")
