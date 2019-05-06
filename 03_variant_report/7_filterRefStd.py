import sys
import pandas as pd



output="/home/hhadmin/exome_pipeline/02_variantCalling/RefStd/"

refstd="/home/hhadmin/exome_pipeline/02_variantCalling/RefStd/somatic_reference_mod.txt"

refstd="/home/hhadmin/environments/ngs_test/exomeAnalysis/Paired/COLO-829_OPT_S5_COLO-829BL_OPT_S4/varscan/COLO-829_OPT_S5_COLO-829BL_OPT_S4.varscan.filter.vcf.txt"


df_ref = pd.read_csv(refstd,sep='\t')

design= "/home/hhadmin/exome_pipeline/agilentCre/cre_design.bed"
# design= "/home/hhadmin/exome_pipeline/01_bamQC/cre_design_exonOnly.bed.txt"
df_design = pd.read_csv(design,sep='\t')
df_design.columns=['CHROM','START','END']

# filter=[]
# for indx,row in df_ref.iterrows():
#     if ( df_design[ (df_design['CHROM'] == row['CHROM']) & (df_design['START'] <= row['POS']) & (df_design['END'] >= row['POS']) ].shape[0]>=1):
#         # print([row['CHROM'],row['POS'],row['REF'],row['ALT'],row['Illumina Alleles'],row['Pleasance Alleles'],row['TGen Alleles'],row['GSC Alleles']])
#         filter.append([row['CHROM'],row['POS'],row['REF'],row['ALT'],row['Illumina Alleles'],row['Pleasance Alleles'],row['TGen Alleles'],row['GSC Alleles']\
#         ,row['Illumina Depth'],row['Pleasance Depth'],row['TGen Depth'],row['GSC Depth']])




filter=[]
for indx,row in df_ref.iterrows():
    if ( df_design[ (df_design['CHROM'] == row['CHROM']) & (df_design['START'] <= row['POS']) & (df_design['END'] >= row['POS']) ].shape[0]>=1):
        print([row['CHROM'],row['POS'],row['REF'],row['ALT']])
        filter.append([row['CHROM'],row['POS'],row['REF'],row['ALT']])

df_filter=pd.DataFrame(filter)
df_filter.columns=['CHROM','POS','REF','ALT','IlluminaAlleles','PleasanceAlleles','TGenAlleles','GSCAlleles'\
,'IlluminaDepth','PleasanceDepth','TGenDepth','GSCDepth']
df_filter.to_csv(output+"reference_standard_creDesign_mod.txt",sep='\t',index=False)



### explore allele frequency and depth of somatic reference
refstd="/home/hhadmin/exome_pipeline/02_variantCalling/RefStd/reference_standard_creDesign_mod.txt"
df_ref = pd.read_csv(refstd,sep='\t')


filter_depth=50
filt_index = df_ref[ (  ( df_ref['IlluminaDepth']<filter_depth) | \
                        ( df_ref['PleasanceDepth']<filter_depth) |\
                        ( df_ref['TGenDepth']<filter_depth) | \
                        ( df_ref['GSCDepth']<filter_depth) ) ].index

df_ref.drop(filt_index,inplace=True)
df_ref.to_csv(output+"reference_standard_creDesign_mod_depth_50.txt",sep='\t',index=False)


def checkAllele(val):
    vals= [int(x)for x in val.split(',')]
    normal=vals[0]
    tumor=vals[1]

    if normal==0:
        return 1
    else :
        n_freq=float(normal)/float(tumor)
        v_freq=float(tumor)/float(normal)
        if n_freq <= 0.5 :#and v_freq >=0.1:
        # if v_freq >=0.1:
            return 1
        else:
            return 0

filter=[]

for indx,row in df_ref.iterrows():

    temp=[]
    temp.append(checkAllele(row['IlluminaAlleles']))
    temp.append(checkAllele(row['PleasanceAlleles']))
    temp.append(checkAllele(row['TGenAlleles']))
    temp.append(checkAllele(row['GSCAlleles']))

    filter.append(temp)


df_ref['Illumina_allele_chk'] = [x[0] for x in filter]
df_ref['Pleasance_allele_chk'] = [x[1] for x in filter]
df_ref['TGen_allele_chk'] = [x[2] for x in filter]
df_ref['GSC_allele_chk'] = [x[3] for x in filter]


filt_index = df_ref[ ( ( df_ref['Illumina_allele_chk'] == 0) | ( df_ref['Pleasance_allele_chk']== 0) | ( df_ref['TGen_allele_chk']== 0) | ( df_ref['GSC_allele_chk']== 0) ) ].index
df_ref.drop(filt_index,inplace=True)
df_ref.to_csv(output+"reference_standard_creDesign_mod_allele.txt",sep='\t',index=False)
