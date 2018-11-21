from fpdf import FPDF
import pandas as pd
import sys

def add_image(pdf, newpage,image_path,title):
    if newpage:
        pdf.add_page()
    pdf.ln(5)  # move 5 down
    pdf.write(5,title)
    pdf.ln(5)
    pdf.image(image_path,w=200,h=100)

def add_text(pdf,newpage,file_path,title):
    if newpage:
        pdf.add_page()
    pdf.ln(5)  # move 5 down
    pdf.write(5,title)
    pdf.ln(10)
    file = open(file_path,'r')
    for line in file:
        pdf.write(5,line)
    file.close()


sample = sys.argv[1]

pdf = FPDF()
pdf.set_font("Arial", size=12)

## add sample Information
pdf.add_page()
pdf.ln(5)  # move 5 down
pdf.write(5,"Sample: COV-1")
pdf.ln(5)

## alignment
df_align = pd.read_csv(sample + ".sorted.bam.alignmentMetrics.txt",skiprows=6,sep='\t')
df_align = df_align.iloc[:,[0,1,2,5,6,7]]
for indx,row in df_align.iterrows():
    pdf.write(5,"\n")
    pdf.write(5,str(row).replace('Name: '+str(indx)+', dtype: object',''))
    pdf.write(5,"\n\n")

df_dups = pd.read_csv(sample + ".sorted.rmdups.bam.metrics.txt",skiprows=6,sep='\t',dtype=str)
df_dups['TOTAL_DUPLICATES'] = str(int(df_dups['UNPAIRED_READ_DUPLICATES'])+(int(df_dups['READ_PAIR_DUPLICATES'])*2)+(int(df_dups['READ_PAIR_OPTICAL_DUPLICATES'])*2))
df_dups = df_dups[['UNPAIRED_READ_DUPLICATES','READ_PAIR_DUPLICATES','READ_PAIR_OPTICAL_DUPLICATES','TOTAL_DUPLICATES','ESTIMATED_LIBRARY_SIZE','PERCENT_DUPLICATION']]
for indx,row in df_dups.iterrows():
    pdf.write(5,str(row).replace('Name: 0, dtype: object',''))

pdf.write(5,"\n\n")
df_target = pd.read_csv("output_hs_metrics.txt",skiprows=6,sep='\t',dtype=str)
df_target = df_target[['BAIT_SET','GENOME_SIZE','TARGET_TERRITORY','TOTAL_READS','PCT_SELECTED_BASES','PCT_USABLE_BASES_ON_TARGET',\
       'FOLD_ENRICHMENT', 'ZERO_CVG_TARGETS_PCT', 'FOLD_80_BASE_PENALTY', \
       'PCT_TARGET_BASES_2X', 'PCT_TARGET_BASES_10X', 'PCT_TARGET_BASES_20X',\
       'PCT_TARGET_BASES_30X', 'PCT_TARGET_BASES_40X', 'PCT_TARGET_BASES_50X',\
       'PCT_TARGET_BASES_100X']]
df_target.columns=['TARGET','GENOME_SIZE','TARGET_TERRITORY','TOTAL_READS','PERCENT_BASES_ON_NEAR_TARGET','PERCENT_BASES_ON_TARGET',\
       'FOLD_ENRICHMENT', 'PERCENT_ZERO_CVG_TARGETS', 'FOLD_80_BASE_PENALTY', \
       'PERCENT_TARGET_BASES_2X', 'PERCENT_TARGET_BASES_10X', 'PERCENT_TARGET_BASES_20X',\
       'PERCENT_TARGET_BASES_30X', 'PERCENT_TARGET_BASES_40X', 'PERCENT_TARGET_BASES_50X',\
       'PERCENT_TARGET_BASES_100X']
for indx,row in df_target.iterrows():
    pdf.write(5,str(row).replace('Name: 0, dtype: object',''))

# pdf.output(sample + "_Alignment.pdf")

add_image(pdf,1,sample + "_500x_coverage_QC.png","Target Region Coverage")
add_image(pdf,0,sample + "_X_coverage_QC.png","Target Region Coverage of Interest")
add_image(pdf,0,"Target_coverage_distribution.png","")
add_image(pdf,1,sample + ".sorted.rmdups.filter.bam.idxstats.plot.png","Distribution of mapped reads among chromosomes")
add_image(pdf,1,sample + ".sorted.rmdups.filter.bammapq.1percent.png","Distribution of map quality <1%")
add_image(pdf,0,sample + ".sorted.rmdups.filter.bammapq.100percent.png","Distribution of map quality 100 %")
add_image(pdf,1,sample + ".sorted.rmdups.filter.bamflag.1percent.png","Distribution of fastq flag <1%")
add_image(pdf,0,sample + ".sorted.rmdups.filter.bamflag.100percent.png","Distribution of fastq flag 100 %")
pdf.output(sample + "_Alignment.pdf")
