import pandas as pd 
import scipy.stats as stats
import common_snps_analysis


def assignRegion(df):
    df["region"]=""
    for indx,row in df.iterrows():
        if indx>=266 and indx<=21555:
            df.ix[indx,"region"] = "orf1ab"
        elif indx>=21563 and indx<=25384:
            df.ix[indx,"region"] = "S"
        elif indx>=25393 and indx<=26220:
            df.ix[indx,"region"] = "ORF3a"
        elif indx>=26245 and indx<=26472:
            df.ix[indx,"region"] = "E"
        elif indx>=26523 and indx<=27191:
            df.ix[indx,"region"] = "M"
        elif indx>=27202 and indx<=27387:
            df.ix[indx,"region"] = "ORF6"
        elif indx>=27394 and indx<=27759:
            df.ix[indx,"region"] = "ORF7a"
        elif indx>=27894 and indx<=28259:
            df.ix[indx,"region"] = "ORF8"
        elif indx>=28274 and indx<=29533:
            df.ix[indx,"region"] = "N"
        elif indx>=29558 and indx<=29674:
            df.ix[indx,"region"] = "ORF10"
    return df


def allStrains():

    df = common_snps_analysis.combineGenome()

    df["strain"] = [x.replace("-0","-") for x in df["strain"].values]
    remove_index = df[df["strain"].isin(["MCoV-1255","MCoV-1343","MCoV-743",
    "MCoV-1219","MCoV-12783","MCoV-12765","MCoV-12805","MCoV-12123","MCoV-12291","MCoV-12788","MCoV-12824","MCoV-12867","MCoV-11867"])].index
    df.drop(remove_index,inplace=True)
    df.rename(columns={'strain':'Strain'},inplace=True)

    ref_dir="/home/tmhsxs240/COVID_19/reference/"
    dflog = pd.read_excel(ref_dir+"4_1_Curated_MCOV_MRN_Strains_withDeceasedCol.xlsx")


    reference = "".join(x for x in df[df["Strain"]=="MN908947"].values[0][1:])


    dfjoin = pd.merge(df,dflog[["Strain","DECEASED"]],on="Strain",how="left",indicator=True)
    dfjoin = dfjoin[dfjoin["_merge"]=="both"]

    fe_test = []
    for position in dfjoin.columns.drop(['Strain','DECEASED','_merge']):

        ctab = pd.crosstab(dfjoin[position],dfjoin["DECEASED"])
        ctab = ctab[ctab.index.isin(['a','t','g','c'])]
        
        if ctab.shape[0]==1:
            continue
        
        pos_result = []

        ref = reference[ctab.index.name-266]
        ref_nd = ctab[ctab.index==ref]['N'].values[0]
        ref_d = ctab[ctab.index==ref]['Y'].values[0]

        pos_result.append(position)
        pos_result.append(ref.upper())
        pos_result.append(ref_nd)
        pos_result.append(ref_d)

        ###remove reference
        ctab = ctab[ctab.index!=ref]

        ### sort to get variants with high number first
        ctab = ctab.sort_values('N',ascending=False)
        
        alt_list = list(ctab.index)
        pos_result.append(len(alt_list))

        for alt in alt_list:
            
            alt_nd = ctab[ctab.index==alt]['N'].values[0]
            alt_d = ctab[ctab.index==alt]['Y'].values[0]

            pos_result.append(alt.upper())
            pos_result.append(alt_nd)
            pos_result.append(alt_d)

            oddsratio, pvalue = stats.fisher_exact([[alt_d,alt_nd],[ref_d,ref_nd]])
            pos_result.append(pvalue)
            pos_result.append(oddsratio)
            
        fe_test.append(pos_result)

    df_fetest = pd.DataFrame(fe_test,columns=['position','ref','ref-ND','ref-D','alt_count',
    'var1','var1-ND','var1-D','var1-pval','var1-oddsr',
    'var2','var2-ND','var2-D','var2-pval','var2-oddsr',
    'var3','var3-ND','var3-D','var3-pval','var3-oddsr',
    ])

    df_fetest.set_index("position",inplace=True)
    df_fetest = assignRegion(df_fetest)
    df_fetest.to_csv("genomewide_fexact_results_all_Strains.csv",index=True)

def finduniqueNumber(dflog,dfjoin,position,variant,status):
    selected_strains = dfjoin[( (dfjoin[position]==variant) & (dfjoin["DECEASED"]==status))]["Strain"].values
    return len(dflog[dflog["Strain"].isin(selected_strains)]["MRN"].unique())

def UniquePatients():
    
    df = common_snps_analysis.combineGenome()

    df["strain"] = [x.replace("-0","-") for x in df["strain"].values]
    remove_index = df[df["strain"].isin(["MCoV-1255","MCoV-1343","MCoV-743",
    "MCoV-1219","MCoV-12783","MCoV-12765","MCoV-12805","MCoV-12123","MCoV-12291","MCoV-12788","MCoV-12824","MCoV-12867","MCoV-11867"])].index
    df.drop(remove_index,inplace=True)
    df.rename(columns={'strain':'Strain'},inplace=True)

    ref_dir="/home/tmhsxs240/COVID_19/reference/"
    dflog = pd.read_excel(ref_dir+"4_1_Curated_MCOV_MRN_Strains_withDeceasedCol.xlsx")


    reference = "".join(x for x in df[df["Strain"]=="MN908947"].values[0][1:])


    dfjoin = pd.merge(df,dflog[["Strain","DECEASED"]],on="Strain",how="left",indicator=True)
    dfjoin = dfjoin[dfjoin["_merge"]=="both"]

    fe_test = []
    for position in dfjoin.columns.drop(['Strain','DECEASED','_merge']):

        ctab = pd.crosstab(dfjoin[position],dfjoin["DECEASED"])
        ctab = ctab[ctab.index.isin(['a','t','g','c'])]
        
        if ctab.shape[0]==1:
            continue
        
        pos_result = []

        ref = reference[ctab.index.name-266]
        ref_nd = ctab[ctab.index==ref]['N'].values[0]
        ref_d = ctab[ctab.index==ref]['Y'].values[0]

        pos_result.append(position)
        pos_result.append(ref.upper())
        pos_result.append(ref_nd)
        pos_result.append(ref_d)

        ###remove reference
        ctab = ctab[ctab.index!=ref]

        ### sort to get variants with high number first
        ctab = ctab.sort_values('N',ascending=False)
        
        alt_list = list(ctab.index)
        pos_result.append(len(alt_list))

        for alt in alt_list:
            
            alt_nd = ctab[ctab.index==alt]['N'].values[0]
            alt_d = ctab[ctab.index==alt]['Y'].values[0]

            pos_result.append(alt.upper())
            pos_result.append(alt_nd)
            pos_result.append(alt_d)
            oddsratio, pvalue = stats.fisher_exact([[alt_d,alt_nd],[ref_d,ref_nd]])
            pos_result.append(pvalue)
            pos_result.append(oddsratio)

            alt_nd_unq = alt_d_unq = 0

            if alt_nd > 4 and alt_d > 4 :
                alt_nd_unq = finduniqueNumber(dflog,dfjoin,position,alt,"N")
                alt_d_unq = finduniqueNumber(dflog,dfjoin,position,alt,"Y")

            pos_result.append(alt_nd_unq)
            pos_result.append(alt_d_unq)
            
        fe_test.append(pos_result)

    df_fetest = pd.DataFrame(fe_test,columns=['position','ref','ref-ND','ref-D','alt_count',
    'var1','var1-ND','var1-D','var1-pval','var1-oddsr','var1-ND-U','var1-D-U',
    'var2','var2-ND','var2-D','var2-pval','var2-oddsr','var2-ND-U','var2-D-U',
    'var3','var3-ND','var3-D','var3-pval','var3-oddsr','var3-ND-U','var3-D-U'
    ])

    df_fetest.set_index("position",inplace=True)
    df_fetest = assignRegion(df_fetest)
    df_fetest.to_csv("genomewide_fexact_results_all_strains_uniqueMRN.csv",index=True)


def UniquePatientsCurated():
    
    df = common_snps_analysis.combineGenome()

    df["strain"] = [x.replace("-0","-") for x in df["strain"].values]
    remove_index = df[df["strain"].isin(["MCoV-1255","MCoV-1343","MCoV-743",
    "MCoV-1219","MCoV-12783","MCoV-12765","MCoV-12805","MCoV-12123","MCoV-12291","MCoV-12788","MCoV-12824","MCoV-12867","MCoV-11867"])].index
    df.drop(remove_index,inplace=True)
    df.rename(columns={'strain':'Strain'},inplace=True)

    ref_dir="/home/tmhsxs240/COVID_19/reference/"
    dflog = pd.read_excel(ref_dir+"4_1_Curated_MCOV_MRN_Strains_withDeceasedCol.xlsx")


    reference = "".join(x for x in df[df["Strain"]=="MN908947"].values[0][1:])


    dfjoin = pd.merge(df,dflog[["Strain","DECEASED"]],on="Strain",how="left",indicator=True)
    dfjoin = dfjoin[dfjoin["_merge"]=="both"]

    fe_test = []

    for position in dfjoin.columns.drop(['Strain','DECEASED','_merge']):

        ctab = pd.crosstab(dfjoin[position],dfjoin["DECEASED"])
        ctab = ctab[ctab.index.isin(['a','t','g','c'])]
        
        if ctab.shape[0]==1:
            continue
        

        ref = reference[ctab.index.name-266]
        ref_nd = finduniqueNumber(dflog,dfjoin,position,ref,"N")
        ref_d = finduniqueNumber(dflog,dfjoin,position,ref,"Y")


        ###remove reference
        ctab = ctab[ctab.index!=ref]

        ### sort to get variants with high number first
        ctab = ctab.sort_values('N',ascending=False)
        
        alt_list = list(ctab.index)
        

        for alt in alt_list:
            
            alt_nd = ctab[ctab.index==alt]['N'].values[0]
            alt_d = ctab[ctab.index==alt]['Y'].values[0]

            if alt_nd > 4 and alt_d > 4 :

                alt_nd_unq = finduniqueNumber(dflog,dfjoin,position,alt,"N")
                alt_d_unq = finduniqueNumber(dflog,dfjoin,position,alt,"Y")

                if alt_nd_unq > 4 and alt_d_unq > 4 :

                    pos_result = []
                    pos_result.append(position)
                    pos_result.append(ref.upper())
                    pos_result.append(ref_nd)
                    pos_result.append(ref_d)

                    oddsratio, pvalue = stats.fisher_exact([[alt_d_unq,alt_nd_unq],[ref_d,ref_nd]])

                    pos_result.append(alt.upper())
                    pos_result.append(alt_nd_unq)
                    pos_result.append(alt_d_unq)

                    pos_result.append(pvalue)
                    pos_result.append(oddsratio)
        
                    fe_test.append(pos_result)

    df_fetest = pd.DataFrame(fe_test,columns=['position','ref','ref-ND','ref-D',
    'var','var-ND','var-D','pval','oddsr'
    ])

    df_fetest.set_index("position",inplace=True)
    df_fetest = assignRegion(df_fetest)
    df_fetest.to_csv("genomewide_fexact_results_all_strains_uniqueMRN_curated.csv",index=True)

# UniquePatients()
UniquePatientsCurated()

## strains per MRN
## 1471 MRNs have 2 strains
# 1     12967
# 2      1471
# 3       330
# 4       108
# 5        52
# 6        18
# 7         8
# 8         8
# 9         6
# 10        4
# 11        2
