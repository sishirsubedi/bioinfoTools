import pandas as pd
from Bio import AlignIO, SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import numpy as np
import os


donor_rec_pair = {"MCoV-35":["MCoV-137","MCoV-230","MCoV-398","MCoV-350","MCoV-449","MCoV-781","MCoV-785"],
            "MCoV-29":["MCoV-101","MCoV-292","MCoV-234","MCoV-349","MCoV-498"],
            "MCoV-11":["MCoV-348","MCoV-347"],
            "MCoV-05":["MCoV-272"]}


alignment_file= "data/4_15/S.fa"
protein_name="S"
protein_reference="MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFR\
SSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIR\
GWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVY\
SSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQ\
GFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFL\
LKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITN\
LCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCF\
TNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYN\
YLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPY\
RVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFG\
RDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAI\
HADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPR\
RARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTM\
YICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFG\
GFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFN\
GLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQN\
VLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGA\
ISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMS\
ECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAH\
FPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELD\
SFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELG\
KYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSE\
PVLKGVKLHYT"

filt = []

ref = []
ref.append("Reference_"+protein_name)
for aa in protein_reference:
    ref.append(aa)

filt.append(ref)

alignment = AlignIO.read(alignment_file,"fasta")

for line in alignment:
    temp =[]
    temp.append(line.id)
    for nt in str(line.seq):
        temp.append(nt)
    filt.append(temp)

df = pd.DataFrame(filt)

selected = []
selected.append("Reference_"+protein_name)
for key, values in donor_rec_pair.items():
    selected.append(key)
    for v in values:
        selected.append(v)

df = df[df[0].isin(selected)]


### write excel
writer = pd.ExcelWriter('excel_file_name.xlsx')
for key, values in donor_rec_pair.items():
     df_subset = df[df[0]=="Reference_"+protein_name]
     df_subset = df_subset.append(df[df[0]==key])
     for v in values:
         df_subset = df_subset.append(df[df[0]==v])
     df_subset.to_excel(writer, sheet_name=key,index=False)
writer.save()
writer.close()

clade={
"Reference_S": "_",
"MCoV-35":	"A2a",
"MCoV-137":	"A2a",
"MCoV-230":	"A2a",
"MCoV-398":	"A2a",
"MCoV-350":	"A2a",
"MCoV-449":	"A2a",
"MCoV-29":	"B1",
"MCoV-101":	"A2a",
"MCoV-292":	"A2a",
"MCoV-234":	"B",
"MCoV-349":	"A2a",
"MCoV-11":	"A2a",
"MCoV-348":	"A2a",
"MCoV-347":	"A2a",
"MCoV-05":	"B1",
"MCoV-272":	"A2a"
}

#### write fasta
for key, values in donor_rec_pair.items():
     df_subset = df[df[0]=="Reference_"+protein_name]
     df_subset = df_subset.append(df[df[0]==key])
     for v in values:
         df_subset = df_subset.append(df[df[0]==v])

     filt = []
     for indx,row in df_subset.iterrows():
         info = row[0]+"-"+clade[row[0]]
         temp =SeqRecord(Seq("".join([x for x in row[1:].values])),info,info,info)
         filt.append(temp)
     align1 = MultipleSeqAlignment(filt)
     AlignIO.write(align1, key+"_alignment.fa","fasta")
