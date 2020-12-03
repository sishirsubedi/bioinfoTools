import pandas as pd
import sys

new_run = sys.argv[1]#"8_"
region_name = sys.argv[2] #S or nsp12
mode = sys.argv[3]# nta or vcf
new_column = sys.argv[4]

# new_run = "8_11"
# region_name = "S"
# mode = "nta"

##############################################################################

out_dir = "data/"+new_run+"/results/"+region_name+"/"
OUT_FILE = out_dir+'COVID_SNP_MAP_'+region_name+'_'+mode+'_table_'+new_run+'.xlsx'
OUT_FILE_all = out_dir+'COVID_SNP_MAP_'+region_name+'_'+mode+'_table_'+new_run+'_all.xlsx'

WORD_FILE=out_dir+region_name+'_table.xlsx'
SNP_FILE=out_dir+'COVID_SNP_MAP_'+region_name+'_'+mode+'_'+new_run+'.xlsx'

# new_column = 'Total Confirmed (1153)'
##########################################################################


df = pd.read_excel(WORD_FILE)
df.columns = ['Genomic locus', 'Gene locus', 'AA Change', 'Domain','Total Confirmed (5085)']
df['Genomic locus'] = df['Genomic locus'].astype(int)
df['Gene locus'] = [x.replace(" ","") for x in df['Gene locus']]

df_snp = pd.read_excel(SNP_FILE)

if mode=="vcf":
    df_snp.columns = ['Genomic locus', 'REF', 'ALT','ALT_count',
            'Annotation', 'Annotation_Impact', 'HGVS.c', 'Gene locus','HGVS.p','AA Change' ,
            'ALL_Variants','Total_Mismatch_Count','Count_per_Variant', 'Variants_Info' ]
elif mode=="nta":
    df_snp.columns = ['Genomic locus', 'REF', 'ALT','ALT_count',
            'Annotation', 'Annotation_Impact', 'HGVS.c', 'Gene locus','HGVS.p','AA Change' ,
            'ALL_Variants','Total_Mismatch_Count', 'Variants_Info' ]

df_snp['Genomic locus'] = df_snp['Genomic locus'].astype(int)
df_snp['Gene locus'] = [x.replace(" ","") for x in df_snp['Gene locus']]

dfjoin = pd.merge(df,df_snp, how='outer',on=['Genomic locus','Gene locus'],indicator=True)
dfjoin[new_column] = dfjoin["ALT_count"]
dfjoin["Variant status"] = ""
dfjoin.ix[dfjoin["_merge"]=="left_only",new_column] = ""
dfjoin.ix[dfjoin["_merge"]=="right_only","AA Change_x"] = dfjoin.ix[dfjoin["_merge"]=="right_only","AA Change_y"].values
dfjoin.ix[dfjoin["_merge"]=="left_only","Variant status"] = "past only"
dfjoin.ix[dfjoin["_merge"]=="right_only","Variant status"] = "current only"
dfjoin.ix[dfjoin["_merge"]=="both","Variant status"] = "past and current"


#### change to syn
# for indx,row in dfjoin[dfjoin["Variant status"]=="current only"].iterrows():
#     aa_change = row['AA Change_x']
#     if aa_change[0] == aa_change[len(aa_change)-1]:
#         dfjoin.ix[indx,'AA Change_x'] = "Syn"


for indx,row in dfjoin.iterrows():
    aa_change = row['AA Change_x']
    if aa_change =="Syn":
        continue
    if aa_change[0] == aa_change[len(aa_change)-1]:
        dfjoin.ix[indx,'AA Change_x'] = "Syn"

protein_domain = {}

if region_name =="S":
    protein_domain={
    'S1 - NTD':'16-305',
    'S1 - RBD':'330-521',
    'Furin Cleavage':'682-685',
    'S2 - FP':'816-833',
    'S2 - HR1':'908-985',
    'S2 - CH':'986-1035',
    'S2 - CD':'1076-1141'}
elif region_name =="nsp12":
    protein_domain={
    'N-terminus':'0-28',
    'B hairpin':'29-50',
    'NiRAN':'51-249',
    'Interface':'250-365',
    'Fingers':'366-581',
    'Palm':'582-620',
    'Fingers':'621-679',
    'Palm':'680-815',
    'Thumb':'816-932'
    }

### add domain
for indx,row in dfjoin.iterrows():
    protein_locus = row['AA Change_x']
    if protein_locus == "Syn":
        continue
    else:
        loc = int(protein_locus[1:len(protein_locus)-1])


        if region_name =="S":
            if loc<=681:
                dfjoin.ix[indx,'Domain'] = 'S1'
            elif loc>=686:
                dfjoin.ix[indx,'Domain'] = 'S2'

        for d in protein_domain:
            start = int(protein_domain[d].split('-')[0])
            end = int(protein_domain[d].split('-')[1])

            if loc>=start and loc<=end:
                dfjoin.ix[indx,'Domain'] = d


dfjoin['AA Change'] = dfjoin["AA Change_x"]
new_columns = [x for x in dfjoin.columns[0:7]]
new_columns[2]="AA Change"
new_columns.append(new_column)
new_columns.append("Variant status")
dfjoin = dfjoin[new_columns]
dfjoin.sort_values(by=['Genomic locus'],inplace=True)

dfjoin.drop(['REF','ALT'], axis=1, inplace=True)

# dfjoin.to_excel(,index=False)

def writeExcel(df,OUT_FILE):
    pd.formats.format.header_style = None
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(OUT_FILE, engine='xlsxwriter')

    # Convert the dataframe to an XlsxWriter Excel object.
    df.to_excel(writer, sheet_name='Sheet1',index=False)

    # Get the xlsxwriter workbook and worksheet objects.
    workbook  = writer.book
    worksheet = writer.sheets['Sheet1']

    font_fmt = workbook.add_format({'align': 'center','font_name': 'Arial', 'font_size': 16,'bold':False})
    header_fmt = workbook.add_format({'align': 'center','text_wrap': True,'font_name': 'Arial', 'font_size': 14, 'bold': True})

    worksheet.set_column('A:G', None, font_fmt)
    worksheet.set_row(0, None, header_fmt)
    worksheet.set_column('A:G', 25)

    number_rows = len(df.index) + 1

    format1 = workbook.add_format({'bg_color': '#FFFF00'})

    worksheet.conditional_format("$A$1:$I$%d" % (number_rows),
                                 {"type": "formula",
                                  "criteria": '=INDIRECT("G"&ROW())="current only"',
                                  "format": format1
                                 }

    )

    writer.save()


writeExcel(dfjoin,OUT_FILE)
