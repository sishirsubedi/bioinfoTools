import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import pandas as pd
import seaborn as sns
from Bio import AlignIO, SeqIO, Entrez
import numpy as np
import os

class Mutation:
    def __init__(self, pos, ref, alt):
        self.pos = pos
        self.ref = ref
        self.alt = alt

class CladeState:
    def __init__(self, name, parent, child,mutations):
        self.name = name
        self.parent = parent
        self.child = child
        self.mutations = mutations

    def addChild(self, new_child):
        self.child = new_child

    def addParent(self, new_parent):
        self.parent = new_parent

class Clade:
    def __init__(self, name):
        self.name = name
        self.clades =[]

    def createNextStrainClade(self):

        clade_19A_mutation_1 = Mutation("8782","","C")
        clade_19A_mutation_2 = Mutation("14408","","C")
        clade_19A = CladeState("19A",None,None,[clade_19A_mutation_1,clade_19A_mutation_2])

        clade_19B_mutation_1 = Mutation("8782","","T")
        clade_19B_mutation_2 = Mutation("28144","","C")
        clade_19B = CladeState("19B",None,None,[clade_19B_mutation_1,clade_19B_mutation_2])

        # clade_20A_mutation_1 = Mutation("8782","C","C")
        clade_20A_mutation_2 = Mutation("14408","","T")
        clade_20A_mutation_3 = Mutation("23403","","G")
        clade_20A = CladeState("20A",None,None,[clade_20A_mutation_2,clade_20A_mutation_3])

        # clade_20B_mutation_1 = Mutation("8782","C","C")
        clade_20B_mutation_2 = Mutation("14408","","T")
        clade_20B_mutation_3 = Mutation("23403","","G")
        clade_20B_mutation_4 = Mutation("28881","","A")
        clade_20B_mutation_5 = Mutation("28882","","A")
        clade_20B = CladeState("20B",None,None,[clade_20B_mutation_2,clade_20B_mutation_3,clade_20B_mutation_4,clade_20B_mutation_5])

        # clade_20C_mutation_1 = Mutation("8782","C","C")
        clade_20C_mutation_2 = Mutation("14408","","T")
        clade_20C_mutation_3 = Mutation("23403","","G")
        clade_20C_mutation_4 = Mutation("25563","","T")
        clade_20C_mutation_5 = Mutation("1059","","T")
        clade_20C = CladeState("20C",None,None,[clade_20C_mutation_2,clade_20C_mutation_3,clade_20C_mutation_4,clade_20C_mutation_5])

        self.clades.append(clade_19A)
        self.clades.append(clade_19B)
        self.clades.append(clade_20A)
        self.clades.append(clade_20B)
        self.clades.append(clade_20C)

        return self.clades

    def createGISAIDClade(self):

        clade_19A_mutation_1 = Mutation("8782","","C")
        clade_19A_mutation_2 = Mutation("28144","","T")
        clade_19A = CladeState("L",None,None,[clade_19A_mutation_1,clade_19A_mutation_2])

        clade_19B_mutation_1 = Mutation("8782","","T")
        clade_19B_mutation_2 = Mutation("28144","","C")
        clade_19B = CladeState("S",None,None,[clade_19B_mutation_1,clade_19B_mutation_2])

        clade_19B_mutation_1 = Mutation("11083","","T")
        clade_19B_mutation_2 = Mutation("26144","","T")
        clade_19B = CladeState("V",None,None,[clade_19B_mutation_1,clade_19B_mutation_2])

        clade_20A_mutation_1 = Mutation("241","","T")
        clade_20A_mutation_2 = Mutation("3037","","T")
        clade_20A_mutation_3 = Mutation("23403","","G")
        clade_20A = CladeState("G",None,None,[clade_20A_mutation_2,clade_20A_mutation_3])

        clade_20B_mutation_1 = Mutation("241","","T")
        clade_20B_mutation_2 = Mutation("3037","","T")
        clade_20B_mutation_3 = Mutation("23403","","G")
        clade_20B_mutation_4 = Mutation("25563","","T")
        clade_20B = CladeState("GH",None,None,[clade_20B_mutation_1,clade_20B_mutation_2,clade_20B_mutation_3,clade_20B_mutation_4])

        clade_20C_mutation_1 = Mutation("241","","T")
        clade_20C_mutation_2 = Mutation("3037","","T")
        clade_20C_mutation_3 = Mutation("23403","","G")
        clade_20C_mutation_4 = Mutation("28882","","A")
        clade_20C = CladeState("GR",None,None,[clade_20C_mutation_1,clade_20C_mutation_2,clade_20C_mutation_3,clade_20C_mutation_4])

        self.clades.append(clade_19A)
        self.clades.append(clade_19B)
        self.clades.append(clade_20A)
        self.clades.append(clade_20B)
        self.clades.append(clade_20C)

        return self.clades

def assignClades(vcfs_directory):

    # vcfs = []
    # for vcfs_directory in vcfs_directory_list:
    #     new_vcfs = [vcfs_directory + x for x in os.listdir(vcfs_directory)]
    #     vcfs = vcfs + new_vcfs


    ### depend on where the files are

    currentClades = Clade("currentClades")
    allclades = currentClades.createGISAIDClade()

    vcfs = os.listdir(vcfs_directory)

    df_combine = pd.DataFrame()
    failed = 0
    for vcf_file in vcfs:
        try:
            df = pd.read_csv(vcfs_directory+vcf_file,sep="\t",skiprows=9)
            df["STRAIN"] = vcf_file.split(".")[0]
            df["Primary_Call"] = [x.split(";")[0].split("=")[1] for x in df.INFO.values]
            df = df[["STRAIN","POS","REF","ALT","QUAL","Primary_Call"]]
            df_combine = pd.concat([df_combine,df],axis=0)
        except:
            failed += 1
            print("failed number:"+str(failed)+"---"+vcf_file)


    ## remove repeat variants
    df_combine["STRAIN"] = [x.replace("-r1","").replace("_r1","").replace("r1","") for x in df_combine["STRAIN"].values]
    df_combine["STRAIN"] = [x.replace("-r2","").replace("_r2","").replace("r2","") for x in df_combine["STRAIN"].values]
    df_combine["STRAIN"] = [x.replace("-r3","").replace("_r3","").replace("r3","") for x in df_combine["STRAIN"].values]
    df_combine["STRAIN"] = [x.replace("v3","").replace("V3","") for x in df_combine["STRAIN"].values]
    df_combine["STRAIN"] = [x.replace("_","-") for x in df_combine["STRAIN"].values]
    df_combine["STRAIN"] = [x.replace("-0","-") for x in df_combine["STRAIN"].values]

    df_combine.drop_duplicates(["STRAIN","POS","REF","ALT"],inplace=True)
    print("Total strains removing repeats : " + str(len(df_combine["STRAIN"].unique())))


    #### filters
    df_combine = df_combine[df_combine["ALT"] == df_combine["Primary_Call"]]

    #filter quality
    df_combine["QUAL"] = df_combine["QUAL"].astype(float)
    df_combine = df_combine[df_combine["QUAL"]>=50.0]

    strain_clades = []
    for strain in df_combine["STRAIN"].unique():
        df_strain = df_combine[df_combine["STRAIN"]==strain]
        found_clades = []
        for clade in allclades:
            found_mutations = 0
            need_mutations = len(clade.mutations)
            for mutation in clade.mutations:
                # if df_strain[( (df_strain["POS"]==int(mutation.pos)) & (df_strain["REF"]==mutation.ref) & (df_strain["ALT"]==mutation.alt) )].shape[0]==1:
                if df_strain[( (df_strain["POS"]==int(mutation.pos)) & (df_strain["ALT"]==mutation.alt) )].shape[0]==1:
                    found_mutations += 1
            if found_mutations == need_mutations:
                found_clades.append(clade.name)
        if len(found_clades)==2 and "G" in found_clades and "GH" in found_clades:
            found_clades = ["GH"]
        elif len(found_clades)==2 and "G" in found_clades and "GR" in found_clades:
            found_clades = ["GR"]
        strain_clades.append([strain,found_clades])
        print(strain,found_clades)


    df_res = pd.DataFrame(strain_clades)
    df_res.columns=["Strain","Clade"]
    df_res.to_excel("strain_clade_assignment.xlsx",index=False)



assignClades("/net/fs01.cluster.com/home/tmhsxs240/COVID_19/data/Ns_Strain_Comparison/6_15_vcfs/")
