import pandas as pd
import os
from pathlib import Path
from more_itertools import chunked
from gen_utils.gen_io import bashCommunicator


### get consensus file and generate alignment in batches
## input all FASTA consensus files
## batch alignment 
def generate_batch_alignment(REFERENCE,file_dir,out_dir):
    
    all_aligned_fasta_files = []
    all_fasta_files = []
    for fasta_files in Path(file_dir).glob('*.fasta'):
        all_fasta_files.append(fasta_files.name)

    df_all_strains = pd.DataFrame(all_fasta_files)

    CHUNK_SIZE = 500

    index_chunks = chunked(df_all_strains.index, CHUNK_SIZE)

    counter = 0
    for ii in index_chunks:
        counter += 1
        for fasta_file in df_all_strains.iloc[ii].values:
            cmd = "cat %s >> %s " % (file_dir+fasta_file[0],out_dir+str(counter)+"_merged.fasta")
            bashCommunicator(cmd)

        cmd = "cat %s  %s > %s " % (out_dir+REFERENCE+".fasta",
                                    out_dir+str(counter)+"_merged.fasta",
                                    out_dir+str(counter)+"_MN908947_merged.fasta")
        
        bashCommunicator(cmd)

        cmd = " mafft --thread 24 --reorder %s > %s " % (out_dir+str(counter)+"_MN908947_merged.fasta",
                                                        out_dir+str(counter)+"_MN908947_merged.mafft_algn.fasta")

        print("running mafft aligner ..."+str(counter))
        bashCommunicator(cmd)

        all_aligned_fasta_files.append(str(counter)+"_MN908947_merged.mafft_algn.fasta")
    
    return all_aligned_fasta_files
