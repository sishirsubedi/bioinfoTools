
import pandas as pd
# import os
import sys

# clinvar = open(sys.argv[1], 'r')
# f = open("clinvar_mod.vcf", "w")
#
# for line in clinvar:
# 	if line[0] == '#':
# 		f.write(line)
# 	else:
# 		elements = line.rstrip().split('\t')
# 		chr="chr"+elements[0]
# 		pos=elements[1]
# 		id=elements[2]
# 		#ref=elements[3]
# 		#alt=elements[4]
# 		info= elements[7]
# 		gene=""
# 		exon=""
# 		info = info.split(';')
# 		for word in info:
# 			if word[0:4]=="GENE":
# 				geneline=word[9:].split(':')
# 				gene=geneline[0]
# 			if word[0:2]=="MC":
# 				wordline=word.split('|')
# 				exon=wordline[1]
# 		if gene=="":gene="nan"
# 		if exon =="": exon="nan"
# 		f.write('%s\t%s\t%s\t%s\t%s\n' %(chr,pos,id,gene,exon))
# f.close()

#df_clinvar = pd.read_csv(sys.argv[1], sep='\t', header=0, skiprows=27)
df_clinvar = pd.read_csv("clinvar_mod.vcf", sep='\t', header=0, skiprows=27)
print(len(df_clinvar.index))
filter_indx = df_clinvar[df_clinvar['CHROM'].str.contains("MT|NW_00")].index
print (len(filter_indx))
df_clinvar.drop(filter_indx, inplace=True)
print(len(df_clinvar.index))
print(df_clinvar.tail())


#df_bedfile = pd.read_csv(sys.argv[2], sep='\t', header=None, skiprows=2)
df_bedfile = pd.read_csv("v7_S31285117_Regions.bed", sep='\t', header=None, skiprows=2)
df_bedfile.columns =['v7_CHROM','v7_START','v7_STOP']
#print (df_bedfile.columns)
#print(df_bedfile.head())

#df_bedfile2 = pd.read_csv(sys.argv[3], sep='\t', header=None, skiprows=2)
df_bedfile2 = pd.read_csv("cre2_S30409818_Regions.bed", sep='\t', header=None, skiprows=2)
df_bedfile2 = df_bedfile2.iloc[:,0:3]
df_bedfile2.columns =['cre2_CHROM','cre2_START','cre2_STOP']
#print (df_bedfile2.columns)
#print(df_bedfile2.head())


genes=df_clinvar['GENE'].unique()
counter =0
combine=[]

for gene in genes:

	counter += 1
	print (counter)

	if gene=="nan":continue

	df_gene = df_clinvar[df_clinvar['GENE']==gene]
	if(len(df_gene.index)==0):continue

	clinvar_positions = df_gene['POS'].values
	chromosome= df_gene['CHROM'].values[0]

	v7=0
	v7_exon=0
	cre2=0
	cre2_exon=0

	### not including exon
	# for clinvar_position in clinvar_positions:
	# 	match = df_bedfile[ (df_bedfile['v7_START'] <= clinvar_position) & (df_bedfile['v7_STOP'] >= clinvar_position) & (df_bedfile['v7_CHROM'] == chromosome)]
	# 	if (len(match.index)>=1):
	# 		v7 = 1
	# 		break
	#
	#
	# for clinvar_position in clinvar_positions:
	# 	match2 = df_bedfile2[ (df_bedfile2['cre2_START'] <= clinvar_position) & (df_bedfile2['cre2_STOP'] >= clinvar_position) & (df_bedfile2['cre2_CHROM'] == chromosome)]
	# 	if (len(match2.index)>=1):
	# 		cre2 = 1
	# 		break


    ## include exons
	for clinvar_position in clinvar_positions:

		match = df_bedfile[ (df_bedfile['v7_START'] <= clinvar_position) & (df_bedfile['v7_STOP'] >= clinvar_position) & (df_bedfile['v7_CHROM'] == chromosome)]
		if (len(match.index)>=1):
			v7 = 1
			v7_exon += 1

		match2 = df_bedfile2[ (df_bedfile2['cre2_START'] <= clinvar_position) & (df_bedfile2['cre2_STOP'] >= clinvar_position) & (df_bedfile2['cre2_CHROM'] == chromosome)]
		if (len(match2.index)>=1):
			cre2 = 1
			cre2_exon += 1

	combine.append([chromosome,gene,len(clinvar_positions),v7,v7_exon, cre, cre2_exon])


df_combine = pd.DataFrame(combine)
df_combine.columns = ['chromosome', 'gene', 'clinvar_exons','v7','v7_exon', 'cre2', 'cre2_exon']
df_combine.to_csv("clinvar_temp.csv")
