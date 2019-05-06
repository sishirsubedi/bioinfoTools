
import pandas as pd
# import os
import sys

cosmic = open("CosmicCodingMuts80.vcf", 'r')

f = open("cosmic_mod.vcf", "w")

for line in cosmic:
	if line[0] == '#':
		f.write(line)
	else:
		elements = line.rstrip().split('\t')
		chr="chr"+elements[0]
		id=elements[2]
		info= elements[7]
		gene="nan"
		info = info.split(';')
		for word in info:
			if word[0:4]=="GENE":
				geneline=word.split('=')
				gene=geneline[1]
				if '_' in gene:
					geneline2 = gene.split('_')
					gene=geneline2[0]
		f.write('%s\t%s\t%s\n' %(chr,id,gene))

f.close()



df_cosmic = pd.read_csv("cosmic_mod.vcf", sep='\t', header=0, skiprows=13)
print(len(df_cosmic.index))
filter_indx = df_cosmic[df_cosmic['CHROM'].str.contains("MT|NW_00")].index
print (len(filter_indx))
df_cosmic.drop(filter_indx, inplace=True)
print(len(df_cosmic.index))
print(df_cosmic.tail())

df_combine= pd.read_csv("clinvar_comparison2_combine.csv",sep=',')
df_combine.shape
df_combine.head()

df_cosmic_genes = df_cosmic.groupby('GENE').size().rename('count').reset_index()
df_cosmic_genes.shape
df_cosmic_genes.head()
df_cosmic_genes.columns=['gene','cosmic_mutations']


df_combine_cosmic = pd.merge(df_combine, df_cosmic_genes, on='gene',how='left',indicator=True)
df_combine_cosmic.head()
df_combine_cosmic.shape

df_combine_cosmic.to_csv("clinvar_comparison2_combine_cosmic.csv",index=False)
