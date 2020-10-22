import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import pandas as pd
import seaborn as sns
import numpy as np
from multiprocessing import Pool,Manager
########################################################

def pre_processes(df_db):
    df_db = df_db[df_db["RESULT"].isin(["Detected","Presumptive Positive"])]
    df_db = df_db[['MRN','ORDER_ID','COLLECTION_DT']]
    df_db.drop_duplicates("ORDER_ID",inplace=True)
    df_mrn_count = df_db.groupby("MRN").size().rename("count").reset_index()
    mrn_of_interest = df_mrn_count[df_mrn_count["count"]>1]["MRN"].values
    del(df_mrn_count)
    return df_db[df_db["MRN"].isin(mrn_of_interest)]


def days_diff(mrn):
    print(mrn)
    collection_dates = df_db[df_db["MRN"]==mrn]["COLLECTION_DT"].values
    out_queue.put([mrn,(np.max(collection_dates) - np.min(collection_dates)).astype('timedelta64[D]')])

def complex_fuction(mrn): ## not working just for NOTE
    print(mrn)
    #### following will not work because Manager proxy objects are unable to 
    # propagate changes made to (unmanaged) mutable objects inside a container.
    # ns.df.ix[ns.df["MRN"]=="109487158","ORDER_ID"]="CHANGED"
    
    ##### this will also not work bc changes made inside process are not global
        df_db.ix[df_db["MRN"]==mrn,"ORDER_ID"]="CHANGED"

    ### fix is to partition dataframe and concatenate at the end

        def work(data_part):
            change_mrn = ['108978149','039566799']
            data_part.ix[data_part["MRN"].isin(change_mrn),"ORDER_ID"]="CHANGED"
            return data_part

        data_split = np.array_split(df_db, 4)
        pool = Pool(4)
        data = pd.concat(pool.map(work, data_split))
        pool.close()
        pool.join()

    ### conclusion
    # 1. if simple function and df does not need to change then just use global df
    # 2. if complex function and df does not need to change then use manager.namespace for proxy
    # 3. if simple or complex function but df does NEED to change then use df split technique 

if __name__ == '__main__':
 
    ref_dir="/home/tmhsxs240/COVID_19/reference/"
    df_db = pd.read_excel(ref_dir+"1_covid_patient_orders.xlsx")
    df_db2 = pd.read_excel(ref_dir+"2_Surveillance.xlsx")
    df_db = df_db.append(df_db2)

    df_db = pre_processes(df_db)

    

    mgr = Manager()

    # we can pass common dataframe using manager namespace, but it makes copy and sends to each process
    # which is not fast for the simple function, it may be faster for other complex situation
    ns = mgr.Namespace()
    ns.df = df_db    


    ## about queue -
    '''
    multiprocessing.Queue() is an object whereas multiprocessing.Manager().Queue() is an address (proxy)
    pointing to shared queue managed by the multiprocessing.Manager() object.
    therefore you can't pass normal multiprocessing.Queue() objects to Pool methods, because it can't be pickled.
    '''
    out_queue=mgr.Queue()

    pool = Pool(processes=4)

    # pool.map(days_diff, df_db.MRN.unique())
    pool.map(complex_fuction, df_db.MRN.unique())

    pool.close()
    pool.join()
    
    results =[]
    for i in range(out_queue.qsize()):
		temp=[]
		temp=out_queue.get()
		results.append(temp)
