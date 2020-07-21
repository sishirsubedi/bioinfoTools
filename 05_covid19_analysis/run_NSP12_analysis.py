import pandas as pd
import SNPanalyzer

old_run = "4_15"
new_run = "7_15"

region_name = "nsp12"
out_dir = "results/"+new_run+"/"+region_name+"/"

reported = "results/"+old_run+"/"


####
vcfs_directory = "data/"+new_run+"/"+new_run+"_vcfs/"
fasta_directory = "data/"+new_run+"/"+new_run+"_fastas/"
SNPanalyzer.analyzeFastaVariants(fasta_directory,vcfs_directory,out_dir)
##############################


vcfs_directory = "data/"+new_run+"/"+new_run+"_vcfs/"
vcf_region=[13442,16236]
SNPanalyzer.analyzeVCF(vcfs_directory,vcf_region,region_name,out_dir,reported)

#
#
# dna_alignment_file = "data/"+new_run+"/Houston.July1.clean--RedundantMRN.fa"
# dna_alignment_region=[13387,16182]
# SNPanalyzer.analyzeGenomicRegion(dna_alignment_file,dna_alignment_region,region_name,out_dir,reported)


# protein_alignment_file = "data/"+new_run+"/RNA-dependent_RNA_polymerase.fa"
# protein_name="nsp12"
# protein_reference="SADAQSFLNVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKTNCCRFQEKDEDDNLIDSYFVVKRHTFSNYQHEETIYNLLKDCPAV\
# AKHDFFKFRIDGDMVPHISRQRLTKYTMADLVYALRHFDEGNCDTLKEILVTYNCCDDDYFNKKDWYDFVENPDILRVYANLGERVRQALLKTVQFCDAMRNAGIVGVLTLD\
# NQDLNGNWYDFGDFIQTTPGSGVPVVDSYYSLLMPILTLTRALTAESHVDTDLTKPYIKWDLLKYDFTEERLKLFDRYFKYWDQTYHPNCVNCLDDRCILHCANFNVLFSTV\
# FPPTSFGPLVRKIFVDGVPFVVSTGYHFRELGVVHNQDVNLHSSRLSFKELLVYAADPAMHAASGNLLLDKRTTCFSVAALTNNVAFQTVKPGNFNKDFYDFAVSKGFFKEG\
# SSVELKHFFFAQDGNAAISDYDYYRYNLPTMCDIRQLLFVVEVVDKYFDCYDGGCINANQVIVNNLDKSAGFPFNKWGKARLYYDSMSYEDQDALFAYTKRNVIPTITQMNL\
# KYAISAKNRARTVAGVSICSTMTNRQFHQKLLKSIAATRGATVVIGTSKFYGGWHNMLKTVYSDVENPHLMGWDYPKCDRAMPNMLRIMASLVLARKHTTCCSLSHRFYRLA\
# NECAQVLSEMVMCGGSLYVKPGGTSSGDATTAYANSVFNICQAVTANVNALLSTDGNKIADKYVRNLQHRLYECLYRNRDVDTDFVNEFYAYLRKHFSMMILSDDAVVCFNS\
# TYASQGLVASIKNFKSVLYYQNNVFMSEAKCWTETDLTKGPHEFCSQHTMLVKQGDDYVYLPYPDPSRILGAGCFVDDIVKTDGTLMIERFVSLAIDAYPLTKHPNQEYADV\
# FHLYLQYIRKLHDELTGHMLDMYSVMLTNDNTSRYWEPEFYEAMYTPHTVLQ"
# SNPanalyzer.analyzeProtein(protein_alignment_file,protein_name,protein_reference,out_dir)

# df_nt = pd.read_excel(result_dir+"COVID_SNP_MAP_"+region_name+"_nt_variants.xlsx")
# df_vcf = pd.read_excel(result_dir+"COVID_SNP_MAP_"+region_name+"_vcfs.xlsx")
# df_aa = pd.read_excel(result_dir+"COVID_SNP_MAP_"+region_name+"_aa_variants.xlsx")
# df_combine = pd.merge(df_nt, df_vcf, how='outer', on=["index"], suffixes=('_nt', '_vcf'))
# df_combine = pd.merge(df_combine, df_aa, how='outer', left_on=["gene position mod"], right_on=["index"], suffixes=('_nt', '_aa'))
# df_combine.to_excel(result_dir+protein_name+"_combine_variants.xlsx",index=False)



# import pandas as pd
# df = pd.read_excel("S/COVID_SNP_MAP_S_NTA_variants.xlsx")
# df.head()
# df.NTA_Status.value_counts()
# df = df[df.NTA_Status=="NotReported"]
# df.shape
# nsp12 = pd.read_csv("sprotein.csv")
# nsp12.head()
# df[df.NTA_Genomic_Locus.isin(nsp12["Genomic locus"].values)]["NTA_Genomic_Locus"].shape
# df[~df.NTA_Genomic_Locus.isin(nsp12["Genomic locus"].values)]["NTA_Genomic_Locus"].shape
# df[~df.NTA_Genomic_Locus.isin(nsp12["Genomic locus"].values)]["NTA_Genomic_Locus"].values
# %history
