import pandas as pd
from boxsdk import Client, OAuth2

def setup_con():

    # Define client ID, client secret, and developer token.
    # get this information from your account on web - https://app.box.com/developers/console
    # you need to create new app for python with oauth2 
    CLIENT_ID = ""
    CLIENT_SECRET = ""
    ACCESS_TOKEN = ""

    oauth2 = OAuth2(CLIENT_ID, CLIENT_SECRET, access_token=ACCESS_TOKEN)

    client = Client(oauth2)

    user = client.user().get()
    print('The current user ID is {0}'.format(user.id))

    return client

def getbam(client,out_dir):
    shared_folder = client.get_shared_item("")

    df = pd.read_csv("")
    strains = df["isolates"].values
    strains = [x.split('.')[0] for x in strains]

    counter =0
    for item in shared_folder.get_items(limit=1000):
        strain_name = item.name.split(".")[0]
        strain_name = strain_name.replace("r1","").replace("r2","").replace("r3","")
        strain_name = strain_name.replace("v1","").replace("v2","").replace("v3","")
        if strain_name in strains and "primertrimmed." in item.name:
            print(item.name,strain_name)
            counter += 1
            print('{0} {1} is named "{2}"'.format(item.type.capitalize(), item.id, item.name))
            with open(out_dir+strain_name+".primertrimmed.sorted.bam", 'wb') as open_file:
                client.file(item.id).download_to(open_file)
                open_file.close()

def getfasta(client,out_dir):

    shared_folder = client.get_shared_item("")

    df = pd.read_csv("")
    strains = df["isolates"].values
    strains = [x.split('.')[0] for x in strains]

    for item in shared_folder.get_items(limit=1000):

        strain_name = item.name.split(".")[0]

        if strain_name not in strains: 
            with open(out_dir+strain_name+".fasta", 'wb') as open_file:
                client.file(item.id).download_to(open_file)
                open_file.close()
