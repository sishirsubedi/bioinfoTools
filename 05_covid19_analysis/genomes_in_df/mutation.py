import pandas as pd

def calculateMutations(var_list):
    mut_counter = 0
    mut_list = []
    for var in var_list:
        if var[1] != var[2]:
            mut_counter += 1
            mut_list.append(var[1].upper()+str(var[0])+var[2].upper())
    return (mut_counter,mut_list)

def genomeWideMutationTable(df):

    REFERENCE = "MN908947"
    df2 = df.T
    df2.rename(columns=df2.iloc[0],inplace=True)
    df2=df2.iloc[1:,:]

    mismatch =[]
    for ci,column in enumerate(df2.columns[1:]):
        df_sel = df2[[REFERENCE,column]][df2[column].isin(['a','t','g','c'])]
        mut_counter,mut_list = calculateMutations(zip(df_sel.index,df_sel[REFERENCE],df_sel[column]))            
        mismatch.append([column,mut_counter,mut_list])

    df_mutationTable = pd.DataFrame(mismatch)
    df_mutationTable.columns = ["Strain","mutation_count","mutation_info"]
    return df_mutationTable
