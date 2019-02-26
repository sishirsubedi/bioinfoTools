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


# sample = sys.argv[1]

pdf = FPDF()
pdf.set_font("Arial", size=12)

sample1="COV-1_COV-2"
sample2="COV-1R_COV-2R"
## add sample Information
pdf.add_page()
pdf.ln(5)  # move 5 down
pdf.write(5,"Sample: " + sample1)
pdf.ln(5)
pdf.write(5,"Sample: " + sample2)
pdf.ln(5)

sample1dir="/home/environments/ngs_test/exomeAnalysis/COV-1_COV-2/varscan/"
sample2dir="/home/environments/ngs_test/exomeAnalysis/COV-1R_COV-2R/varscan/"


# ## alignment
# df_align = pd.read_csv(sample + ".sorted.bam.alignmentMetrics.txt",skiprows=6,sep='\t')
# df_align = df_align.iloc[:,[0,1,2,5,6,7]]
# for indx,row in df_align.iterrows():
#     pdf.write(5,"\n")
#     pdf.write(5,str(row).replace('Name: '+str(indx)+', dtype: object',''))
#     pdf.write(5,"\n\n")
#
# df_dups = pd.read_csv(sample + ".sorted.rmdups.bam.metrics.txt",skiprows=6,sep='\t',dtype=str,nrows=1)
# df_dups['TOTAL_DUPLICATES'] = str(int(df_dups['UNPAIRED_READ_DUPLICATES'])+(int(df_dups['READ_PAIR_DUPLICATES'])*2)+(int(df_dups['READ_PAIR_OPTICAL_DUPLICATES'])*2))
# df_dups = df_dups[['UNPAIRED_READ_DUPLICATES','READ_PAIR_DUPLICATES','READ_PAIR_OPTICAL_DUPLICATES','TOTAL_DUPLICATES','ESTIMATED_LIBRARY_SIZE','PERCENT_DUPLICATION']]
# for indx,row in df_dups.iterrows():
#     pdf.write(5,str(row).replace('Name: 0, dtype: object',''))
#
# pdf.write(5,"\n\n")
# df_target = pd.read_csv(sample + ".output_hs_metrics.txt",skiprows=6,sep='\t',dtype=str)
# df_target = df_target[['BAIT_SET','GENOME_SIZE','TARGET_TERRITORY','TOTAL_READS','PCT_SELECTED_BASES','PCT_USABLE_BASES_ON_TARGET',\
#        'FOLD_ENRICHMENT', 'ZERO_CVG_TARGETS_PCT', 'FOLD_80_BASE_PENALTY', \
#        'PCT_TARGET_BASES_2X', 'PCT_TARGET_BASES_10X', 'PCT_TARGET_BASES_20X',\
#        'PCT_TARGET_BASES_30X', 'PCT_TARGET_BASES_40X', 'PCT_TARGET_BASES_50X',\
#        'PCT_TARGET_BASES_100X']]
# df_target.columns=['TARGET','GENOME_SIZE','TARGET_TERRITORY','TOTAL_READS','PERCENT_BASES_ON_NEAR_TARGET','PERCENT_BASES_ON_TARGET',\
#        'FOLD_ENRICHMENT', 'PERCENT_ZERO_CVG_TARGETS', 'FOLD_80_BASE_PENALTY', \
#        'PERCENT_TARGET_BASES_2X', 'PERCENT_TARGET_BASES_10X', 'PERCENT_TARGET_BASES_20X',\
#        'PERCENT_TARGET_BASES_30X', 'PERCENT_TARGET_BASES_40X', 'PERCENT_TARGET_BASES_50X',\
#        'PERCENT_TARGET_BASES_100X']
# for indx,row in df_target.iterrows():
#     pdf.write(5,str(row).replace('Name: 0, dtype: object',''))
#
# # pdf.output(sample + "_Alignment.pdf")

names= [  \
          '_DIST_PLOT_normal_depth' ,\
          '_DIST_PLOT_normal_depth_2',\
          '_DIST_PLOT_tumor_depth',\
          '_DIST_PLOT_tumor_depth_2',\
          '_DIST_PLOT_normal_var_freq_percent' ,\
          '_DIST_PLOT_normal_var_freq_percent_2',\
          '_DIST_PLOT_tumor_var_freq_percent',\
          '_DIST_PLOT_normal_reads1_plus',\
          '_DIST_PLOT_normal_reads1_minus',\
          '_DIST_PLOT_tumor_reads1_plus',\
          '_DIST_PLOT_tumor_reads1_plus_2',\
          '_DIST_PLOT_tumor_reads1_minus',\
          '_DIST_PLOT_tumor_reads1_minus_2',\
          '_DIST_PLOT_tumor_reads2_plus',\
          '_DIST_PLOT_tumor_reads2_plus_2',\
          '_DIST_PLOT_tumor_reads2_minus',\
          '_DIST_PLOT_tumor_reads2_minus_2',\
        ]


for n in names:
    add_image(pdf,1,sample1dir + sample1 + n + ".png",sample1)
    add_image(pdf,0,sample2dir + sample2 + n + ".png",sample2)


pdf.output(sample1dir + sample1 + sample2 + "_Report.pdf")
