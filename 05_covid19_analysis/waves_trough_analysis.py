import pandas as pd
import NTA_SNP_analyzer
import alignment_analysis


def protein_analysis():
    region_name="nsp12"    
    region=[]

    if region_name == "S":
        region="21508:25330"
    elif region_name == "nsp12":
        region="13387:16182"

    reference_alignment = "/home/tmhsxs240/COVID_19/data/4_15/Houston.4-14.clean.fa"
    all_strains = NTA_SNP_analyzer.sequenceFromAlignment(reference_alignment,region_name,region,reference=True)
    
    ###########add all other alignments
    alignment_file ="/home/tmhsxs240/COVID_19/data/8_11_all_5085/Houston.Aug12.clean.fa"
    all_strains.extend(NTA_SNP_analyzer.sequenceFromAlignment(alignment_file,region_name,region))
    
    alignment_file ="/home/tmhsxs240/COVID_19/data/11_20/Nov-19.trimmed.clean.fa"
    all_strains.extend(NTA_SNP_analyzer.sequenceFromAlignment(alignment_file,region_name,region,adjustment=16))

    alignment_file ="/home/tmhsxs240/COVID_19/data/12_1/Nov-29.trimmed.clean.fa"
    all_strains.extend(NTA_SNP_analyzer.sequenceFromAlignment(alignment_file,region_name,region,adjustment=26))

    df = pd.DataFrame(all_strains)
    df = NTA_SNP_analyzer.filterMainAlignment(df)


    ref_dir="/home/tmhsxs240/COVID_19/reference/"
    dflog = pd.read_excel(ref_dir+"4_1_Curated_MCOV_MRN_Strains.xlsx")

    ###if need to save a copy
    # dfjoin = pd.merge(df.iloc[:,0:2],dflog,left_on=0,right_on="Strain",how='outer',indicator=True)
    # dfjoin.to_excel(ref_dir+"5_1_Current_NTA_Curated_MCOV_MRN_Strains_FINAL.xlsx",index=False)


    out_dir = "/home/tmhsxs240/COVID_19/data/12_1/results/"+region_name+"/"

    for analysis_group in dflog.AnalysisGroup.unique():

        ### select one group at a time 
        
        current_strains = list(dflog[dflog.AnalysisGroup==analysis_group]["Strain"].values)
        current_strains.append("MN908947")
        df_group = df[df[0].isin(current_strains)]
        
        print(analysis_group)
        print(df_group.head())
        print(df_group.shape)
        
        df_variants = NTA_SNP_analyzer.getVariantsTable(df_group,region_name)
        df_annot = NTA_SNP_analyzer.addAnnotation(df_variants,out_dir+analysis_group)
        NTA_SNP_analyzer.curationAfterAnnotation(df_annot,region_name,out_dir+analysis_group)


def protein_analysis():

    analysis_group="5085"
    region_name="nsp12"


    region=[]
    if region_name == "S":
        region="21508:25330"
    elif region_name == "nsp12":
        region="13387:16182"

    reference_alignment = "/home/tmhsxs240/COVID_19/data/4_15/Houston.4-14.clean.fa"
    all_strains = NTA_SNP_analyzer.sequenceFromAlignment(reference_alignment,region_name,region,reference=True)
    
    ###########add all other alignments
    alignment_file ="/home/tmhsxs240/COVID_19/data/8_11_all_5085/Houston.Aug12.clean.fa"
    all_strains.extend(NTA_SNP_analyzer.sequenceFromAlignment(alignment_file,region_name,region))
    
    # alignment_file ="/home/tmhsxs240/COVID_19/data/11_20/Nov-19.trimmed.clean.fa"
    # all_strains.extend(NTA_SNP_analyzer.sequenceFromAlignment(alignment_file,region_name,region,adjustment=16))

    # alignment_file ="/home/tmhsxs240/COVID_19/data/12_1/Nov-29.trimmed.clean.fa"
    # all_strains.extend(NTA_SNP_analyzer.sequenceFromAlignment(alignment_file,region_name,region,adjustment=26))

    df = pd.DataFrame(all_strains)
    df = NTA_SNP_analyzer.filterMainAlignment(df)



    out_dir = "/home/tmhsxs240/COVID_19/data/12_1/results/"+region_name+"/"

    df_variants = NTA_SNP_analyzer.getVariantsTable(df,region_name)
    df_annot = NTA_SNP_analyzer.addAnnotation(df_variants,out_dir+analysis_group)
    NTA_SNP_analyzer.curationAfterAnnotation(df_annot,region_name,out_dir+analysis_group)


# ref_dir="/home/tmhsxs240/COVID_19/reference/"
# alignment_analysis.tableGenomicMutationsMultiple(ref_dir+"6_1_Current_NTA_Curated_MCOV_MRN_Strains_genomic_table.csv")
