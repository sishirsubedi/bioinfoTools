import pandas as pd
import SNPanalyzer

old_run = "4_15"
new_run = "7_15"

region_name = "S"
out_dir = "results/"+new_run+"/"+region_name+"/"

reported = "results/"+old_run+"/"

############################
vcfs_directory = "data/"+new_run+"/"+new_run+"_vcfs/"
fasta_directory = "data/"+new_run+"/"+new_run+"_fastas/"
SNPanalyzer.analyzeFastaVariants(fasta_directory,vcfs_directory,out_dir)
################# s protein

vcfs_directory = "data/"+new_run+"/"+new_run+"_vcfs/"
vcf_region=[21563,25384]
SNPanalyzer.analyzeVCF(vcfs_directory,vcf_region,region_name,out_dir,reported)


# dna_alignment_file = "data/"+new_run+"/Houston.July1.clean--RedundantMRN.fa"
# dna_alignment_region=[21508,25330]
# SNPanalyzer.analyzeGenomicRegion(dna_alignment_file,dna_alignment_region,region_name,out_dir,reported)



# protein_alignment_file = "data/"+new_run+"/RNA-dependent_RNA_polymerase.fa"
# alignment_file= "data/S.fa"
# protein_name="S"
# protein_reference="MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFR\
# SSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIR\
# GWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVY\
# SSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQ\
# GFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFL\
# LKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITN\
# LCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCF\
# TNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYN\
# YLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPY\
# RVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFG\
# RDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAI\
# HADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPR\
# RARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTM\
# YICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFG\
# GFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFN\
# GLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQN\
# VLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGA\
# ISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMS\
# ECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAH\
# FPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELD\
# SFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELG\
# KYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSE\
# PVLKGVKLHYT"
# SNPanalyzer.analyzeProtein(alignment_file,protein_name,protein_reference)

####### combine all results

# df_nt = pd.read_excel(result_dir+"COVID_SNP_MAP_"+region_name+"_nt_variants.xlsx")
# df_vcf = pd.read_excel(result_dir+"COVID_SNP_MAP_"+region_name+"_vcfs.xlsx")
# df_aa = pd.read_excel(result_dir+"COVID_SNP_MAP_"+region_name+"_aa_variants.xlsx")
# df_combine = pd.merge(df_nt, df_vcf, how='outer', on=["index"], suffixes=('_nt', '_vcf'))
# df_combine = pd.merge(df_combine, df_aa, how='outer', left_on=["gene position mod"], right_on=["index"], suffixes=('_nt', '_aa'))
# df_combine.to_excel(result_dir+protein_name+"_combine_variants.xlsx",index=False)
