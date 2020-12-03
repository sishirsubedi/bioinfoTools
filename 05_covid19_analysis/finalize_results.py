import pandas as pd
import sys

##############################################################################


def processVariants(strain_group,region_name,out_dir):
    
    SNP_FILE = out_dir + strain_group + "COVID_SNP_MAP_"+region_name+"_NTA.xlsx"

    df_snp = pd.read_excel(SNP_FILE)

    df_snp.columns = ['Genomic locus', 'REF', 'ALT','ALT_count',
                'Annotation', 'Annotation_Impact', 'HGVS.c', 'Gene locus','HGVS.p','AA Change' ,
                'ALL_Variants','Total_Mismatch_Count', 'Variants_Info' ]

    df_snp = df_snp[['Genomic locus', 'ALT_count','Gene locus','AA Change' ]]

    df_snp['Genomic locus'] = df_snp['Genomic locus'].astype(int)
    df_snp['Gene locus'] = [x.replace(" ","") for x in df_snp['Gene locus']]


    for indx,row in df_snp.iterrows():
        aa_change = row['AA Change']
        if aa_change =="Syn":
            continue
        if aa_change[0] == aa_change[len(aa_change)-1]:
            df_snp.ix[indx,'AA Change'] = "Syn"

    protein_domain = {}

    if region_name =="S":
        protein_domain={
        'S1-NTD':'16-305',
        'S1-RBD':'330-521',
        'Furin Cleavage':'682-685',
        'S2-FP':'816-833',
        'S2-HR1':'908-985',
        'S2-CH':'986-1035',
        'S2-CD':'1076-1141'}
    elif region_name =="nsp12":
        protein_domain={
        'N-terminus':'0-28',
        'B hairpin':'29-50',
        'NiRAN':'51-249',
        'Interface':'250-365',
        'Fingers-N':'366-581',
        'Palm-N':'582-620',
        'Fingers-C':'621-679',
        'Palm-C':'680-815',
        'Thumb':'816-932'
        }

    ### add domain
    for indx,row in df_snp.iterrows():
        protein_locus = row['AA Change']
        if protein_locus == "Syn":
            continue
        else:
            loc = int(protein_locus[1:len(protein_locus)-1])


            if region_name =="S":
                if loc<=681:
                    df_snp.ix[indx,'Domain'] = 'S1'
                elif loc>=686:
                    df_snp.ix[indx,'Domain'] = 'S2'

            for d in protein_domain:
                start = int(protein_domain[d].split('-')[0])
                end = int(protein_domain[d].split('-')[1])

                if loc>=start and loc<=end:
                    df_snp.ix[indx,'Domain'] = d

    df_snp.sort_values(by=['Genomic locus'],inplace=True)
    df_snp = df_snp[['Genomic locus','Gene locus','AA Change','Domain','ALT_count' ]]
    df_snp.rename(columns={"ALT_count":strain_group},inplace=True)
    return df_snp


protein="S"

out_dir = "/home/tmhsxs240/COVID_19/data/12_1/results/"+protein+"/"

dfw1 = processVariants("wave_1",protein,out_dir)

dfw2 = processVariants("wave_2",protein,out_dir)

dftr = processVariants("trough",protein,out_dir)

dfjoin = pd.merge(dfw1,dfw2,on=['Genomic locus', 'Gene locus', 'AA Change', 'Domain'],how='outer',indicator=True).rename(columns={"_merge":"wave1-wave2"})

dfjoin = pd.merge(dfjoin,dftr,on=['Genomic locus', 'Gene locus', 'AA Change', 'Domain'],how='outer',indicator=True).rename(columns={"_merge":"wave-trough"})

dfjoin.rename(columns={"wave_1":"wave_1(1026)","wave_2":"wave_2(6310)","trough":"trough(3933)"},inplace=True)

dfjoin.to_excel(out_dir+"wave_1_2_trough_final_snp_table_"+protein+"_protein.xlsx",index=False)




protein="S"

out_dir = "/home/tmhsxs240/COVID_19/data/12_1/results/"+protein+"/"

dfw1 = processVariants("5085",protein,out_dir)

dfw2 = processVariants("previous",protein,out_dir)

dftr = processVariants("new_run",protein,out_dir)

dfjoin = pd.merge(dfw1,dfw2,on=['Genomic locus', 'Gene locus', 'AA Change', 'Domain'],how='outer',indicator=True).rename(columns={"_merge":"wave1-wave2"})

dfjoin = pd.merge(dfjoin,dftr,on=['Genomic locus', 'Gene locus', 'AA Change', 'Domain'],how='outer',indicator=True).rename(columns={"_merge":"wave-trough"})

dfjoin.rename(columns={"5085":"published(5085)","previous":"previous(5417)","new_run":"new_run(768)"},inplace=True)

dfjoin.drop("wave1-wave2",axis=1,inplace=True)
dfjoin.rename(columns={"wave-trough":"status"},inplace=True)

dfjoin["status"]= ["prior and current" if x=="both" else "prior only" if x=="left_only" else "current only"
if x=="right_only" else None for x in dfjoin["status"].values]

dfjoin.to_excel(out_dir+"published5085_previous_new_run_snp_table_"+protein+"_protein.xlsx",index=False)
