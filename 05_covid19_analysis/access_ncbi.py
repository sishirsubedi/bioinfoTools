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
df.to_csv("ncbi.csv")
