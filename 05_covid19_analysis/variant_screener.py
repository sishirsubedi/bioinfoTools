import pandas as pd
import os
import sys
sys.path.append('')

from alignment.generate_alignment import generate_batch_alignment
from alignment.alignment_tools import get_df_from_multiple_alignments
from genomes_in_df.mutation import genomeWideMutationTable
import nt_to_aa_conversion
from gen_utils.gen_io import generateRefFasta,writesummary
import gen_utils.gen_covid as covid

def checkdeletion(variant_del,strain_del):
    return len([x for x in variant_del if x in strain_del])

def variant_count_mutations(df,variants):
    variant_counter=[]
    for indx,row in df.iterrows():
        variant_counter.append(len([x for x in row["mutation_info"] if x in variants]))
    return variant_counter

def get_variants_match_mutations(df,variants):
    all_variants=[]
    for indx,row in df.iterrows():
        sel_variants = [x for x in row["mutation_info"] if x in variants.keys()]
        if len(sel_variants) >0:
            all_variants.append([ variants[x] for x in sel_variants])
        else:
            all_variants.append([])
    return [len(x) for x in all_variants ],[",".join(x) for x in all_variants]

def b117_count_deletions(df):
    deletion_region_start= 11000
    deletion_region_stop= 22000
    df2 = df.T
    df2.rename(columns=df2.iloc[0],inplace=True)
    df2=df2.iloc[1:,:]
    b117_deletions=[]
    for column in df2.columns[1:]:
        b117_deletions.append(list(df2[df2[column]=='-'].loc[deletion_region_start:deletion_region_stop].index))
    return b117_deletions

def b117_mutation_occurrence_isolates(df,b117_variants,b117_deletions):

    var_dict = {'A23063T': 0,
    'A28111G': 0,
    'A28281T': 0,
    'C23271A': 0,
    'C23604A': 0,
    'C23709T': 0,
    'C27972T': 0,
    'C28977T': 0,
    'C3267T': 0,
    'C5388A': 0,
    'G24914C': 0,
    'G28048T': 0,
    'G28280C': 0,
    'T24506G': 0,
    'T28282A': 0,
    'T6954C': 0,
    'b117_del_11288':0,
    'b117_del_21765':0,
    'b117_del_21991':0
    }

    for indx,row in df.iterrows():
        for variant in [x for x in row["mutation_info"] if x in b117_variants]:
            var_dict[variant] += 1

    for indx,row in df.iterrows():
        for b117_del in b117_deletions:
            if checkdeletion(b117_deletions[b117_del],row["B117_deletions"]) == len(b117_deletions[b117_del]):
                var_dict[b117_del] += 1 
    
    return var_dict

def b117_screen_haplotypes(df,b117_variants,b117_deletions):

    mut_dict = {
    'C3267T':'ORF1ab-T1001I',
    'C5388A':'ORF1ab-A1708D',
    'T6954C': 'ORF1ab-I2230T',
    'A23063T':'spike-N501Y',
    'C23271A':'spike-A570D',
    'C23604A':'spike-P681H',
    'C23709T':'spike-T716I',
    'T24506G':'spike-S982A',
    'G24914C':'spike-D1118H',
    'C27972T':'Orf8-Q27stop',
    'G28048T':'Orf8-R52I',
    'A28111G':'Orf8-Y73C',
    'G28280C':'N-D3L-80',
    'A28281T':'N-D3L-81',
    'T28282A':'N-D3L-82',
    'C28977T':'N-S235F',
    'b117_del_11288':'ORF1ab-SGF-3675-3677-deletion',
    'b117_del_21765':'spike-HV-69-70-deletion',
    'b117_del_21991':'spike-Y144-deletion',
    '':''
    }

    haplotyes = []

    for indx,row in df.iterrows():

        mutations = []
        
        for variant in [x for x in row["mutation_info"] if x in b117_variants]:
            mutations.append(variant)

        for b117_del in b117_deletions:
            if checkdeletion(b117_deletions[b117_del],row["B117_deletions"]) == len(b117_deletions[b117_del]):
                mutations.append(b117_del) 
            
        haplotyes.append([row["Strain"],mutations])
    
    haplos_dict = {}
    for ht in haplotyes:
        haplo = "-".join(x for x in ht[1])

        if haplo in haplos_dict:
            haplos_dict[haplo] += 1
        else:
            haplos_dict[haplo] = 1
    

    haplos_aa_dict = {}
    for ht in haplos_dict:
        isolates = haplos_dict[ht]

        muts = ""
        for m in ht.split('-'):
            muts += mut_dict[m]+";"

        haplos_aa_dict[muts] = isolates
    

    df_haplo = pd.DataFrame.from_dict(haplos_aa_dict,orient='index')
    df_haplo.reset_index(inplace=True)
    df_haplo.columns = ["haplotype","total_isolates"]
    return df_haplo

def b117_mutation_occurrence_haplotype(df):
    all_muts = {
    'ORF1ab-T1001I':0,
    'ORF1ab-A1708D':0,
    'ORF1ab-I2230T':0,
    'spike-N501Y':0,
    'spike-A570D':0,
    'spike-P681H':0,
    'spike-T716I':0,
    'spike-S982A':0,
    'spike-D1118H':0,
    'Orf8-Q27stop':0,
    'Orf8-R52I':0,
    'Orf8-Y73C':0,
    'N-D3L-80':0,
    'N-D3L-81':0,
    'N-D3L-82':0,
    'N-S235F':0,
    'ORF1ab-SGF-3675-3677-deletion':0,
    'spike-HV-69-70-deletion':0,
    'spike-Y144-deletion':0,
    '':0
    }

    for indx,row in df.iterrows():
        for mutation in row["haplotype"].split(";"):
            all_muts[mutation] += 1

    return all_muts

def haplotype_ranking(haplotype):

    mutations = ['ORF1ab-T1001I',
    'ORF1ab-A1708D',
    'ORF1ab-I2230T',
    'spike-N501Y',
    'spike-A570D',
    'spike-P681H',
    'spike-T716I',
    'spike-S982A',
    'spike-D1118H',
    'Orf8-Q27stop',
    'Orf8-R52I',
    'Orf8-Y73C',
    'N-D3L-80',
    'N-D3L-81',
    'N-D3L-82',
    'N-S235F',
    'ORF1ab-SGF-3675-3677-deletion',
    'spike-HV-69-70-deletion',
    'spike-Y144-deletion']

    code =""
    for m in mutations:
        if m in haplotype:
            code += '1'
        else:
            code += '0'
    
    return code

def run_screener(file_dir,out_dir):

    REFERENCE = covid.covid_basics.reference

    generateRefFasta(REFERENCE,out_dir)

    all_aligned_fasta_files = generate_batch_alignment(REFERENCE,file_dir,out_dir)

    df = get_df_from_multiple_alignments(all_aligned_fasta_files,out_dir)

    df_muttable = genomeWideMutationTable(df)

    df_muttable = nt_to_aa_conversion.add_aa_mutations_column(df_muttable,out_dir)

    b117_variants_dict = covid.covid_variants.b117_variants
    b117_deletions = covid.covid_variants.b117_deletions
    cal20_variants_dict = covid.covid_variants.cal20_variants
    b1351_variants_dict = covid.covid_variants.b1351_variants
    b1128_variants_dict = covid.covid_variants.b1128_variants
    p2_variants_dict = covid.covid_variants.p2_variants
    
    df_muttable["CAL20_variants"] ,df_muttable["CAL20_variants_out_of_6"] =\
    get_variants_match_mutations(df_muttable,cal20_variants_dict)

    df_muttable["B117_variants"],df_muttable["B117_variants_out_of_16"]  =\
    get_variants_match_mutations(df_muttable,b117_variants_dict)
    
    df_muttable["B1351_variants"],df_muttable["B1351_variants_out_of_12"] =\
    get_variants_match_mutations(df_muttable,b1351_variants_dict)

    df_muttable["B1128_variants"],df_muttable["B1128_variants_out_of_15"] =\
    get_variants_match_mutations(df_muttable,b1128_variants_dict)

    df_muttable["P2_variants"],df_muttable["P2_variants_out_of_20"] =\
    get_variants_match_mutations(df_muttable,p2_variants_dict)

    df_muttable["B117_deletions"] = b117_count_deletions(df)

    ## save
    df_muttable.to_excel(out_dir+"variants_screen_snps.xlsx",index=False)

def main():
    out_dir = ""
    file_dir = ""
    run_screener(file_dir,out_dir)

if __name__ == "__main__":
    main()