import pandas as pd
import ast
import pylab
import numpy as np
import os
import re

def clonesfromSNPtable(snp_table):
    df= pd.read_excel(snp_table)

    df = df[['NTA_Genomic_Locus', 'NTA_Referene', 'NTA_Mismatch_Count', 'NTA_Gene','NTA_Variants', 'NTA_Variants_Info']]

    combine =[]
    for indx,row in df.iterrows():
        alt_variant = ast.literal_eval(row["NTA_Variants"])[0]
        variants = ast.literal_eval(row["NTA_Variants_Info"])
        for variant in variants:
            if variant == alt_variant:
                combine.append([row["NTA_Mismatch_Count"],row["NTA_Referene"].upper()+str(row['NTA_Genomic_Locus'])+alt_variant.upper(),variants[variant][1]])
    df_combine = pd.DataFrame(combine)
    df_combine.columns =[ 'NTA_Mismatch_Count','NTA_Variants', 'NTA_Variants_Info']

    df_combine["count_len"] = [len(x) for x in df_combine["NTA_Variants_Info"]]

    df_combine  = df_combine[df_combine["count_len"]==df_combine["NTA_Mismatch_Count"]]


    df_combine.NTA_Mismatch_Count.value_counts()
    df_combine = df_combine[df_combine["NTA_Mismatch_Count"]>1]

    clone_match={}
    for num in df_combine.NTA_Mismatch_Count.unique():

        clone_match[num]=0
        df_selected = df_combine[df_combine["NTA_Mismatch_Count"]==num]

        all_strains = []
        for indx,strain_list in df_selected.NTA_Variants_Info.items():
            for strain in strain_list:
                if strain not in all_strains:
                    all_strains.append(strain)

        repeat_strains = {}
        for s in all_strains: repeat_strains[s] =0

        for indx,strain_list2 in df_selected.NTA_Variants_Info.items():
            for strain2 in all_strains:
                if strain2 in strain_list2:
                    repeat_strains[strain2] += 1

        for s in repeat_strains:
            if repeat_strains[s] >1:
                clone_match[num] += 1

    df_combine.to_excel(out_file,index=False)

def clonesfromSNPtable(clones_from_snp_table,out_file):

    df_clonal = pd.read_excel(clones_from_snp_table)

    zip_codes =[]
    mrn_codes =[]
    for indx,row in df_clonal.iterrows():
        zip = []
        mrn = []
        strain_list = ast.literal_eval(row.NTA_Variants_Info)
        for strain in strain_list:
            print(strain)
            strain = strain.replace("-0","-")
            if dfjoin_selected[dfjoin_selected.Strain==strain].shape[0] != 0:
                zip.append(dfjoin_selected[dfjoin_selected.Strain==strain]["ZIP"].values[0])
                mrn.append(dfjoin_selected[dfjoin_selected.Strain==strain]["MRN"].values[0])
            else:
                zip.append(0)
                mrn.append(0)
        zip_codes.append(zip)
        mrn_codes.append(mrn)

    df_clonal["zip_codes"] = zip_codes

    df_clonal["MRNs"] = mrn_codes

    df_clonal["zip_codes_counts"] = [len(set(x)) for x in df_clonal["zip_codes"]]

    df_clonal["MRN_counts"] = [len(set(x)) for x in df_clonal["MRNs"]]

    df_clonal.to_excel(out_file)


def clonesfromSNPtable(clones_from_snp_table,out_file):

    from matplotlib.backends.backend_pdf import PdfPages
    from uszipcode import SearchEngine

    df_clonal = pd.read_excel(clones_from_snp_table)
    BBox = (-96.3940,-94.5813,30.3551,29.3247)

    #####plot zipcodes

    df_clonal_mini = df_clonal[df_clonal.NTA_Mismatch_Count>=0]

    df_clonal_mini.sort_values("NTA_Mismatch_Count",ascending=False,inplace=True)

    search = SearchEngine(simple_zipcode=True)

    with PdfPages(out_file) as pdf:
        for indx,row in df_clonal_mini.iterrows():

            fig, ax = plt.subplots(figsize = (15,10))
            for zcode in row.zip_codes:

                # zcdb = ZipCodeDatabase()

                if zcode is None or zcode=='nan': continue
                print(zcode)
                try:
                    zipcode = search.by_zipcode(int(zcode))

                    ax.scatter(zipcode.lng, zipcode.lat, zorder=1, c='r', s=50)
                except:
                    print("failed---"+str(zcode))


            ax.set_title(row.NTA_Variants+" , total isolates- "+str(row.NTA_Mismatch_Count)+" , total patient- "+str(row.MRN_counts)+ " ,total zipcodes- "+str(row.zip_codes_counts))
            ax.imshow(ruh_m, zorder=0, extent = BBox, aspect= 'equal')
            pdf.savefig(fig)
