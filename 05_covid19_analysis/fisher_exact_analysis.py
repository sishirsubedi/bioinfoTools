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

def fetest():
    
    df = common_snps_analysis.combineGenome()

    fe_test = []

    for position in dfjoin.columns[2:]:

        ctab = pd.crosstab(dfjoin[position],dfjoin["DECEASED"])
        ctab = ctab[ctab.index.isin(['a','t','g','c'])]
        
        ref = reference[ctab.index.name-266]

        ###remove reference
        ctab = ctab[ctab.index!=ref]
        ### sort to get variants with high number first
        ctab = ctab.sort_values('N',ascending=False)
        
        alt_list = list(ctab.index)
        

        for alt in alt_list:
        
            alt_nd_unq_mrns = finduniqueMRNs(dflog,dfjoin,position,alt,"N")
            alt_d_unq_mrns = finduniqueMRNs(dflog,dfjoin,position,alt,"Y")
            ref_nd_mrns = finduniqueMRNs(dflog,dfjoin,position,ref,"N")
            ref_d_mrns = finduniqueMRNs(dflog,dfjoin,position,ref,"Y")

            alt_nd_unq = len(alt_nd_unq_mrns) 
            alt_d_unq = len(alt_d_unq_mrns) 
            ref_nd = len(ref_nd_mrns) 
            ref_d = len(ref_d_mrns) 


            print(position)

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

            pos_result.append(len(np.intersect1d(alt_d_unq_mrns,ref_d_mrns)))
            pos_result.append(len(np.intersect1d(alt_d_unq_mrns,ref_nd_mrns)))
            pos_result.append(len(np.intersect1d(alt_nd_unq_mrns,ref_d_mrns)))
            pos_result.append(len(np.intersect1d(alt_nd_unq_mrns,ref_nd_mrns)))


            alt_nd_unq = alt_nd_unq - len(np.intersect1d(alt_nd_unq_mrns,ref_nd_mrns))
            alt_d_unq = alt_d_unq  - len(np.intersect1d(alt_d_unq_mrns,ref_d_mrns))
            ref_nd = ref_nd  - len(np.intersect1d(alt_nd_unq_mrns,ref_nd_mrns))
            ref_d = ref_d -   len(np.intersect1d(alt_d_unq_mrns,ref_d_mrns))

            pos_result.append(ref_nd)
            pos_result.append(ref_d)
            oddsratio, pvalue = stats.fisher_exact([[alt_d_unq,alt_nd_unq],[ref_d,ref_nd]])
            pos_result.append(alt_nd_unq)
            pos_result.append(alt_d_unq)

            pos_result.append(pvalue)
            pos_result.append(oddsratio)



            fe_test.append(pos_result)

    df_fetest = pd.DataFrame(fe_test,columns=['position','ref',
    'ref-ND','ref-D',
    'var','var-ND','var-D','pval','oddsr','alt_d_ref_d','alt_d_ref_nd','alt_nd_ref_d','alt_nd_ref_nd',
    'ref-ND_adj','ref-D_adj','var-ND_adj','var-D_adj','pval_adj','oddsr_adj'
    ])

    df_fetest.set_index("position",inplace=True)
    df_fetest = assignRegion(df_fetest)
    df_fetest.to_excel("../genomewide_fexact_results.xlsx",index=True)

