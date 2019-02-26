
import pandas as pd


refstd="/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/reference_standard.txt"
design= "/home/hhadmin/exome_pipeline/01_bamQC/cre_design.bed"
# design= "/home/hhadmin/exome_pipeline/01_bamQC/cre_design_exonOnly.bed.txt"
output="/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/"

df_ref = pd.read_csv(refstd,sep='\t')

df_design = pd.read_csv(design,sep='\t')
df_design.columns=['CHROM','START','END']

filter=[]
for indx,row in df_ref.iterrows():
    if ( df_design[ (df_design['CHROM'] == row['CHROM']) & (df_design['START'] <= row['POS']) & (df_design['END'] >= row['POS']) ].shape[0]>=1):
        print([row['CHROM'],row['POS'],row['REF'],row['ALT']])
        filter.append([row['CHROM'],row['POS'],row['REF'],row['ALT']])

df_filter=pd.DataFrame(filter)
df_filter.columns=['CHROM','POS','REF','ALT']
df_filter.to_csv(output+"reference_standard_creDesign.txt",sep='\t',index=False)
