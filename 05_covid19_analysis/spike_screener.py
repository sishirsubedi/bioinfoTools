import pandas as pd
import os
import sys 
sys.path.append('')
from genomes_in_df.mutation import genomeWideMutationTable
from gen_utils.gen_io import writesummary


def get_all_muts(gene_snps_db):

    df_all_muts = pd.read_excel(gene_snps_db)
    df_all_muts["genomic_info"] = [ x[0]+str(y)+x[len(x)-1] for x,y in zip(df_all_muts["Gene locus"],df_all_muts["Genomic locus"])]
    df_all_muts["protein_info"] = [ x+"("+str(y)+")" for x,y in zip(df_all_muts["AA Change"],df_all_muts["Domain"])]
    spike_mutations = {x:y for x,y in zip(df_all_muts["genomic_info"],df_all_muts["protein_info"])}

    spike_mutations['C22452A']="S297Stop(*)"
    spike_mutations['C24245T']="Q895Stop(*)"
    spike_mutations['C21629T']="Q23Stop(*)"
    spike_mutations['A22145T']="K195Stop(*)"

    spike_mutations['G22487T']="E309Stop(*)"
    spike_mutations['C22982T']="Q474Stop(*)"
    spike_mutations['A25193T']="K1211Stop(*)"
    spike_mutations['A22892T']="K444Stop(*)"
    spike_mutations['C21962T']="Q134Stop(*)"
    spike_mutations['G22403T']="E281Stop(*)"
    spike_mutations['C25163T']="Q1201Stop(*)"

    spike_mutations['G22409T']="G283Stop(*)"

    return spike_mutations


def run_spike_variant_screener(region_alignment,gene_snps_db,out_dir):

    df = pd.read_csv(region_alignment)
    df_muttable = genomeWideMutationTable(df)
    spike_mutations = get_all_muts(gene_snps_db)

    if variant != "gene":
        
        sel_index =[]
        for indx,row in df_muttable.iterrows():
            if variant in row["mutation_info"]:
                sel_index.append(indx)
        
        df_muttable = df_muttable.iloc[sel_index,:]

    haplos_dict = {}
    for indx,row in df_muttable.iterrows():
        
        haplo = "-".join(x for x in row['mutation_info'])

        if haplo in haplos_dict:
            haplos_dict[haplo] += 1
        else:
            haplos_dict[haplo] = 1


    haplos_aa_dict = {}
    for ht in haplos_dict:
        isolates = haplos_dict[ht]

        muts = ""
        for m in ht.split('-'):
            aa_change = spike_mutations[m]
            if "Syn" not in aa_change: 
                muts += spike_mutations[m]+";"

        if muts in haplos_aa_dict:
            haplos_aa_dict[muts] += isolates
        else:
            haplos_aa_dict[muts] = isolates

    writesummary(haplos_aa_dict,out_dir+variant+"_spike-mutation_background.csv")

def main():
    
    region_alignment=""
    gene_snps_db=""
    out_dir = ""
    run_spike_variant_screener(region_alignment,gene_snps_db,out_dir)

if __name__ == "__main__":
    main()


