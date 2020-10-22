import pandas as pd
import SNP_analyzer_main
import sys

old_run = sys.argv[1]#"4_15"
new_run = sys.argv[2]#"8_11"
region_name = sys.argv[3] ## nsp12 or S 
mode = sys.argv[4]
wave = sys.argv[5]

out_dir = "data/"+new_run+"/results/"+region_name+"/"
reported = "results/"+old_run+"/"


print(old_run,new_run,mode)
###########################

def run_vcf_analysis():

    vcfs_directory = "data/"+new_run+"/"+new_run+"_vcfs/"
    fasta_directory = "data/"+new_run+"/"+new_run+"_fastas/"
    SNP_analyzer_main.analyzeFastaVariants(fasta_directory,vcfs_directory,out_dir)

    vcfs_directory = "data/"+new_run+"/"+new_run+"_vcfs/"

    vcf_region = []
    if region_name == "S":
        vcf_region=[21563,25384]
    elif region_name == "nsp12":
        vcf_region=[13442,16236]

    SNP_analyzer_main.analyzeVCF(vcfs_directory,vcf_region,region_name,out_dir,reported)


def run_nta_analysis():
    dna_alignment_file = "data/"+new_run+"/Houston.Sept.clean.fa"

    dna_alignment_region=[]
    if region_name == "S":
        dna_alignment_region=[21508,25330]
    elif region_name == "nsp12":
        dna_alignment_region=[13387,16182]

    SNP_analyzer_main.analyzeGenomicRegion(dna_alignment_file,dna_alignment_region,region_name,out_dir,reported)

def run_nta_analysis_wave(wave):
    dna_alignment_file = "data/"+new_run+"/Houston.Aug12.clean.fa"

    dna_alignment_region=[]
    if region_name == "S":
        dna_alignment_region=[21508,25330]
    elif region_name == "nsp12":
        dna_alignment_region=[13387,16182]

    SNP_analyzer_main.analyzeGenomicRegionWave(dna_alignment_file,dna_alignment_region,region_name,out_dir,reported,wave)


if mode=="vcf":
    run_vcf_analysis()
elif mode=="nta":
    run_nta_analysis()
elif mode=="nta-wave":
    run_nta_analysis_wave(wave)
