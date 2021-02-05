import pandas as pd
from Bio import SeqIO,AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment 
from pathlib import Path
import glob
import NTA_SNP_analyzer

'''

here is one approach using genome assemblies and pairwise alignment with the reference 
( align with nuleotide to find start position of protein,
translate and align with protein to find deletion)

WT: AIHVSG  and GVYYHK 
MUT: AI--SG for HV-69-70 deletion and GV-YHK for Y144 deletion

Found many variations but none for HV-69-70 deletion and 3 strains for Y144 deletion

- download genomes from All Files/..

- for each strain
    1. find start position of Spike protein in genome
       by aligning with the first 30bp position of spike protein reference (nulceotide)
    2. find end position of Spike protein based on the length of spike protein (nulceotide)
    3. extract protein gene sequence based on above start and stop position
    4. translate nuleotide
    5. pairwise alignment of aminoacid sequence with the reference 
    6. check sequences in the alignment at specified positions for HV-69-70 deletion and Y144 deletion
'''

region_name = "S"
region="21508:25330"
reference_alignment = ""
reference_strain = NTA_SNP_analyzer.sequenceFromAlignment(reference_alignment,region_name,region,reference=True)
reference_sequence = "".join(x for x in reference_strain[0][1:])
refseq = Seq(reference_sequence)
refprot = refseq.translate()


final_list =[]
file_dir=""
input_files = Path(file_dir).glob('*.fa')

for input_file in input_files:
    print(input_file.as_posix())
    fasta_sequences = SeqIO.parse(open(input_file.as_posix()),'fasta')

    for strain in fasta_sequences:

        try:
            name, sequence = strain.id.split("/")[2].replace("TX-HMH-",""), str(strain.seq.lower())
            # print(name)

            # if name in ["MCoV-10027","MCoV-11908","MCoV-9551"]:
            # if name in ["MCoV-13133"]:
            start = sequence.find(reference_sequence[0:30])
            stop = start + len(reference_sequence)
            strain_sequence = Seq(sequence[start:stop])
            strainprot = strain_sequence.translate()

            alignments = pairwise2.align.globalxx(refprot, strainprot)

            for a in alignments:
                mut1_pos = str(a).split(',')[0].find("AIHVSG")
                mut2_pos = str(a).split(',')[0].find("GVYYHK")
                final_list.append([name,str(a).split(',')[1][mut1_pos:mut1_pos+6],str(a).split(',')[1][mut2_pos:mut2_pos+6]])
                break
        except :
            print("error")
df = pd.DataFrame(final_list)
df.to_csv("del_test.v2.csv",index=False)
