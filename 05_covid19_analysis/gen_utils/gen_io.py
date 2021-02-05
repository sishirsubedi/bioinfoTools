import os
import sys


def bashCommunicator(command,output_expected=False):

    import subprocess

    process = subprocess.Popen([command],shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print("Process failed %s %d %s %s" % (command,process.returncode, stdout, stderr))
    else:
        if output_expected:
            return [x for x in stdout.split("\n")]


def generateRefFasta(virus,out_dir):

    from Bio import Entrez

    filename = virus+".fasta"
    Entrez.email = "sample@example.com"  
    if not os.path.isfile(filename):
        net_handle = Entrez.efetch(db="nucleotide", id=virus, rettype="fasta", retmode="text")
        out_handle = open(out_dir+filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print("Saved")


def referenceAttributes(out_file):

    import pandas as pd
    from Bio import AlignIO, SeqIO, Entrez

    #### get nsp12 protein sequence from ncbi
    Entrez.email = "sample@example.org"
    handle = Entrez.efetch(db="nucleotide", id="MN908947.3", rettype="gb",retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    gene_location = []
    for feature in record.features:
        if feature.type == "CDS":
            location = feature.location
            gene = feature.qualifiers['gene']
            gene_location.append([gene[0],location.start.position,location.end.position])
    df=pd.DataFrame(gene_location)
    df.to_csv(out_file)


def writesummary(summary,filename):

    import csv

    with open(filename, "w") as f:
        writer = csv.writer(f)
        for i in summary:
            writer.writerow([i, summary[i]])
    f.close() 



