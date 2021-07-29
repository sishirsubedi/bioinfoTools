import pandas as pd
import numpy as np
from datetime import timedelta
import scipy.optimize as optim
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from gen_utils.gen_io import read_run_params,log_msg



#############################################

params = read_run_params()
run = params["current_run"]
out_home = params["container"]+"output/" 
out_dir = out_home+run+"/"

df = pd.read_csv(out_dir+"4_mcov_strain_variant_map_covid_pangolin_db_input_"+run+".csv")
df = df[df.quality=="HQ"]


     
#########################
tag="B.1.617.Family"
voi=["B.1.617.2","AY.2","AY.3"]
start_date = "4-15-2021"
end_date = "7-20-2021"
days_since="4/15/2021"
days= 180

# voi="P.1"
# start_date = "1-1-2021"
# end_date = "6-20-2021"
# days_since="1/1/2021"
# days= 360
#################################


###take unique patients with variant
keep_mrns_variant = np.unique(df[df.variant.isin(voi)]["MRN"])
df_mrns = df[df.MRN.isin(keep_mrns_variant)]
df_mrns = df_mrns[df_mrns.variant.isin(voi)] ###important step--remove non b117 variant 
df_mrns.sort_values("COLLECTION_DT",inplace=True)
df_mrns.drop_duplicates("MRN",keep="first",inplace=True)


keep_mrns_not_variant = np.unique(df[~df.variant.isin(voi)]["MRN"])
df_mrns_not_variant = df[df.MRN.isin(keep_mrns_not_variant)]
df_mrns_not_variant = df_mrns_not_variant[~df_mrns_not_variant.variant.isin(voi)]
df_mrns_not_variant.sort_values("COLLECTION_DT",inplace=True)
df_mrns_not_variant.drop_duplicates("MRN",keep="first",inplace=True)

df_2 = df_mrns.append(df_mrns_not_variant)
df_2.drop_duplicates("MRN",keep="first",inplace=True)

df = df_2


df=df[['MCoVNumber','COLLECTION_DT','variant']]

#####################################

df.COLLECTION_DT = pd.to_datetime(df.COLLECTION_DT)
df.COLLECTION_DT = df.COLLECTION_DT.dt.date


df = df[ (  (df.COLLECTION_DT>=pd.to_datetime(start_date)) &
            (df.COLLECTION_DT<pd.to_datetime(end_date)) 
         )
       ]
df.sort_values("COLLECTION_DT",inplace=True)

df.variant.fillna(0,inplace=True)
#########################

df.variant = [1 if x in voi else 0 for x in df.variant]


df_variant = df.groupby("COLLECTION_DT")["variant"].agg("sum").reset_index()
df_count = df.groupby("COLLECTION_DT")["variant"].agg("count").reset_index()

dates = pd.date_range(df.COLLECTION_DT.min(), (df.COLLECTION_DT.max() + timedelta(days=1) )-timedelta(days=1),freq='d')
df_data = pd.DataFrame(dates)
df_data.columns=["dates"]
df_data["date_step"]= [x for x in range(1,df_data.shape[0]+1,1)]
df_data["total"] = df_count.variant
df_data["variant"] = df_variant.variant
df_data["variant_csum"] = np.cumsum(df_variant.variant.values)
df_data["variant_percent"]=[ (x/y)*100 for x,y in zip(df_data.variant,df_data.total)]
df_data.to_excel("final_Data_"+tag+"_log_growth_6_28_2021.xlsx",index=False)

def my_logistic(x,a,b,c):
     return c/(1 + a * np.exp(-b*x))

x = np.array(df_data.date_step)
# y = np.array(df_data.variant_csum)
y = np.array(df_data.variant_percent)

##########optimize
po = np.random.exponential(size=3)
bounds = (0,[1000.,2.0,100.])
(a,b,c),cov = optim.curve_fit(my_logistic,x,y,bounds=bounds,p0=po)

# for i in range(1,20,1):
#      try:
#           # po = np.array([250.,0.10,99.])
#           po= np.random.exponential(size=3)
#           bounds = ([0.,0.1,0.],[1000.,float(i),100.])
#           (a,b,c),cov = optim.curve_fit(my_logistic,x,y,bounds=bounds,p0=po)
#           print(c)
#      except:
#           print("error for  " + str(i))

# po = np.array([250.,0.10,99.])
# bounds = ([0.,0.1,99.],[1000.,1.0,100.])
# (a,b,c),cov = optim.curve_fit(my_logistic,x,y,bounds=bounds,p0=po)

plt.scatter(x,y)
plt.plot(x,my_logistic(x,a,b,c))
xprime = np.array([x for x in range(1,170,1)])
yprime = my_logistic(xprime,a,b,c)
plt.plot(xprime,yprime)
plt.savefig("log_fit_best_fit"+tag+".png")
plt.close()


############################## method 2 using t distribution on error --> perfer this one 

from scipy.stats.distributions import  t

pars, pcov = (a,b,c),cov

alpha = 0.05 # 95% confidence interval = 100*(1-alpha)

n = len(y)    # number of data points
p = len(pars) # number of parameters

dof = max(0, n - p) # number of degrees of freedom

# student-t value for the dof and confidence level
tval = t.ppf(1.0-alpha/2., dof) 

val_dw = 0
val_up = 0
for i, p,var in zip(range(n), pars, np.diag(pcov)):
    sigma = var**0.5
    
    if i==1:
         val_dw = p - sigma*tval
         val_up = p + sigma*tval

    print ('p{0}: {1} [{2}  {3}]'.format(i, p,
                                  p - sigma*tval,
                                  p + sigma*tval))



plt.plot(x,y,'bo',markersize=5,label='Observed')
xprime = np.array([x for x in range(1,days,1)])
yprime = my_logistic(xprime,a,b,c)
plt.plot(xprime,yprime,label='Predicted')

xpred = np.array([x for x in range(1,days,1)])
ypred_dw = my_logistic(xpred,pars[0],val_dw,pars[2])
ypred_up = my_logistic(xpred,pars[0],val_up,pars[2])

plt.fill_between(xpred, ypred_up,ypred_dw,color = 'k', alpha = 0.1,label='95% CI')

plt.title("Logistic growth model ["+tag+"]",fontsize=18)
plt.xlabel("Days since "+days_since,fontsize=15)
plt.ylabel("Percent of patients ",fontsize=15)

plt.legend()
plt.savefig("log_pred_best_fit"+tag+".png")
plt.close()


gr=b;dt = 70/(gr*100);print(dt)
gr=val_up;dt = 70/(gr*100);print(dt)
gr=val_dw;dt = 70/(gr*100);print(dt)

