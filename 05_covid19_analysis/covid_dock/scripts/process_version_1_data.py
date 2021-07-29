import pandas as pd
from gen_utils.gen_io import read_run_params,log_msg
from alignment.alignment_tools import genomicSequenceFromAlignment
from genomes_in_df.filter import filterMainAlignment

def get_combine_df():

    params = read_run_params()

    all_strains = []

    nta_dir = params["container"]+"input_alignment/"
    
    nta_files = params["nta_group_v1"]

    for nta_group in sorted(nta_files) :
        log_msg("run group is - "+ nta_group)
        align_file = nta_files[nta_group].split(",")[0]
        adj = int(nta_files[nta_group].split(",")[1])

        if "reference" in nta_group:
            all_strains.extend(genomicSequenceFromAlignment(nta_dir+align_file,reference=True))
        else:            
            all_strains.extend(genomicSequenceFromAlignment(nta_dir+align_file,adjustment=adj))

    df = pd.DataFrame(all_strains)
    df = filterMainAlignment(df,"cds")

    return df
