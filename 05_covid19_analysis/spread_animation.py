import pandas as pd
import sys
import plotly.express as px
import json
from multiprocessing import Process
import subprocess

def bashCommunicator(command,output_expected=False):
    process = subprocess.Popen([command],shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print("Process failed %s %d %s %s" % (command,process.returncode, stdout, stderr))
    else:
        print("bash command completed.")
        if output_expected:
            return [x for x in stdout.split("\n")]

def get_zipcode_boundary_geojson(geojson_file):
    ###read city boundary data
    with open(geojson_file) as data:
        geojson = json.load(data)
    return geojson

def get_zipcodes(zipcodes_file):
    selected_counties =['Austin','Brazoria','Chambers','Fort Bend','Galveston','Harris','Liberty','Montgomery','Waller']
    ###get houston zipcodes from https://www.zipcodestogo.com/Texas/
    df_texax_zipcodes = pd.read_excel(zipcodes_file)
    return df_texax_zipcodes[df_texax_zipcodes.County.isin(selected_counties)]["Zip Code"].values

def read_metadata(meta_file,mode):
    df = pd.read_excel(meta_file)
    if mode =="all":
        df = df[df.IS_FIRST_FOR_PATIENT==1]
    else:
        df = df[df[mode] != "NotPresent"]
        #### select unique MRN
        df = df.drop_duplicates('MRN')


    df = df[['COLLECTION_DATE','ZIP_CODE','Musser Lab No.']]
    df = df[df.ZIP_CODE.apply(lambda x: x.isnumeric())]
    df.columns=['COLLECTION_DATE', 'ZCTA5CE10', 'STRAIN']
    df_zipdata = df.groupby(['COLLECTION_DATE','ZCTA5CE10']).count()
    df_zipdata.reset_index(inplace=True)

    ######
    color_max = 50
    if mode != "all":
        color_max = df_zipdata.STRAIN.max()
    ######

    return (df_zipdata,color_max)

def extract_data_for_zipcodes(city_zipcodes,df_zipdata):
    all_data=[]
    all_houston_zipcode_dict = {int(x):0 for x in city_zipcodes}
    for d in df_zipdata.COLLECTION_DATE.unique():
        df_day = df_zipdata[df_zipdata.COLLECTION_DATE==d]
        day_zipcode_dict = {int(x):int(y) for x,y in zip (df_day.ZCTA5CE10.values,df_day.STRAIN.values)}
        for zcode in day_zipcode_dict.keys():
            if zcode in all_houston_zipcode_dict.keys():
                all_houston_zipcode_dict[zcode] += day_zipcode_dict[zcode]
        all_data.append([str(d)[:10],all_houston_zipcode_dict.copy()])
    return all_data

def generate_image(all_data_instance,data_instance,geojson,mode,color_max):
    day = all_data_instance[0]
    zipc_dict = all_data_instance[1]

    df_allzipcode_perday = pd.DataFrame(zipc_dict.items(), columns=['ZCTA5CE10', 'Cumulative Patients'])
    df_allzipcode_perday["Day"] = day
    fig = px.choropleth(df_allzipcode_perday, geojson=geojson, color="Cumulative Patients",
                locations="ZCTA5CE10", featureidkey="properties.ZCTA5CE10",
                projection="mercator",
                range_color=(0, color_max),
                color_continuous_scale=px.colors.sequential.YlOrRd
            )
    ### update_goes to include only locations of interest i.e houston zipcodes
    fig.update_geos(fitbounds="locations", visible=False)
    fig.update_layout(
         title={
        'text': 'Cumulative SARS-CoV-2 strains sequenced <br>'+ mode + " : "+ day,
        'y':0.95,
        'x':0.5,
        'xanchor': 'center',
        'yanchor': 'top'},
         margin={"r":0,"t":0,"l":0,"b":0}
    )

    filename=""
    if len(str(data_instance))==1:
        filename="00"+str(data_instance)
    elif len(str(data_instance))==2:
        filename="0"+str(data_instance)
    elif len(str(data_instance))==3:
        filename=str(data_instance)

    fig.write_image(mode+"/"+ filename+".png",scale=3.0)
    print("Completed- "+filename+"---"+day)

def process_images(all_data,geojson,mode,color_max,p_count):
    for indx in range(0,len(all_data),p_count):
        processes = []
        if indx < len(all_data)-1:
            for i in range(0,p_count):
                process = Process(target=generate_image,args=(all_data[indx+i],indx+i,geojson,mode,color_max))
                processes.append(process)
                process.start()
        else:##if length of data is odd then avoid index error for last element + 1
            generate_image(all_data[len(all_data)-1],len(all_data)-1,geojson,mode,color_max)

        for p in processes:
            p.join()

def create_directory(mode):
    ###create directory
    # cmd = "mkdir /home/covid/covid_paper_animation/%s " %(mode )
    cmd = "mkdir %s " %(mode )

    bashCommunicator(cmd)

def create_video(mode):
    cmd = "ffmpeg -r 5 -f image2 \
    -s 1920x1080 -i /home/covid/covid_paper_animation/%s/%s \
    -vcodec libx264 -crf 25 -pix_fmt yuv420p \
    /home/covid/covid_paper_animation/%s/%s.mp4"%(mode,'"%03d.png"',mode,mode)
    bashCommunicator(cmd)

houston_geojson = get_zipcode_boundary_geojson("../tx_texas_zip_codes_geo.min.json")
houston_zipcodes = get_zipcodes("../texas_zipcodes.xlsx")

######
# df_mutations = pd.read_csv("mutations_input_2.csv")
# df_mutations = pd.read_csv(sys.argv[1])
######

# for mode in df_mutations.mutation:

#     print("processing...."+mode)

#     create_directory(mode)

#     df_data_per_zip,color_max = read_metadata("MCoV Sample Log 8-14-20 metadata_biallelic_column_fixed.xlsx",mode)

#     #### get cumulative data
#     cumulative_zip_data = extract_data_for_zipcodes(houston_zipcodes,df_data_per_zip)

#     print("Number of images to generate...."+str(len(cumulative_zip_data)))

#     process_images(cumulative_zip_data,houston_geojson,mode,color_max,2)

#     create_video(mode)


df = pd.read_excel("postion_mcov_mrn.xlsx")
df.sort_values("COLLECTION-DATE",inplace=True)
df["COLLECTION-DATE"] =[x[0:10] for x in df["COLLECTION-DATE"].values]
print(df.head())

for pos in df.POS.unique():

    if pos in [23756,24076]:

        print(str(pos))

        create_directory(pos)

        
        all_houston_zipcode_dict = {int(x):0 for x in houston_zipcodes}
        df_zipdata = df[df.POS==pos]

        REF=df_zipdata["REF"].unique()[0]
        ALT=df_zipdata["ALT"].unique()[0]

        counter =0
        for day in df_zipdata["COLLECTION-DATE"].unique():

            df_day = df_zipdata[df_zipdata["COLLECTION-DATE"]==day]

            all_data=[]
            for zcode in df_day.ZIP.values:
                zcode = int(str(zcode)[0:5])
                try:
                    all_houston_zipcode_dict[zcode] += 1
                except:
                    print("not found"+str(zcode))
            all_data.append([day,all_houston_zipcode_dict.copy()])

            zipc_dict = all_data[0][1]

            zip_row =[]
            max_val = 0
            for key,val in zipc_dict.items():
                zip_row.append([key,val])
                if max_val < val:
                    max_val = val 

            df_allzipcode_perday = pd.DataFrame(zip_row, columns=['ZCTA5CE10', 'Cumulative Patients'])
            df_allzipcode_perday["Day"] = day
            fig = px.choropleth(df_allzipcode_perday, geojson=houston_geojson, color="Cumulative Patients",
                        locations="ZCTA5CE10", featureidkey="properties.ZCTA5CE10",
                        projection="mercator",
                        range_color=(0, max_val),
                        color_continuous_scale=px.colors.sequential.YlOrRd
                    )
            ### update_goes to include only locations of interest i.e houston zipcodes
            fig.update_geos(fitbounds="locations", visible=False)
            fig.update_layout(
                title={
                'text': 'Cumulative SARS-CoV-2 strains sequenced <br>'+ REF+str(pos)+ALT + " : "+ day,
                'y':0.95,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'},
                margin={"r":0,"t":0,"l":0,"b":0}
            )


            fig.write_image(str(pos)+"/00"+ str(counter)+".png",scale=3.0)
            print("Completed- "+day+"--pic-number--"+str(counter))
            counter +=1


    
