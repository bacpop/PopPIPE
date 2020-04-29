import os, sys
import pandas as pd
from datetime import date
import requests
import json    

# Constants
microreact_api_url = "https://microreact.org/api/project/"
description_string = "PopPIPE run on " + date.today().strftime("%Y-%b-%d %H:%M")

# Format data
with open(snakemake.input['tree'], 'r') as tree_file:
    tree_string = tree_file.read()

with open(snakemake.input['dot'], 'r') as dot_file:
    dot_string = dot_file.read()

clusters = pd.read_table(snakemake.input['clusters'], sep=",",).set_index("Taxon")
clusters.rename(columns={i: i + "__autocolour" for i in clusters.columns})
csv_string = clusters.to_csv(None, sep=",", header=True, index=True)

# Make request
payload = {"name": snakemake.params['microreact_name'],
           "description": description_string,
           "email": snakemake.params['microreact_email'],
           "website": snakemake.params['microreact_website'],
           
           "data": csv_string,
           "tree": tree_string,
           "dot": dot_string} 

r = requests.post(microreact_api_url, data = json.dumps(payload))
if r.status_code == 200:
    with open(snakemake.output, 'w') as url_file:
        url_file.write(r.url)
        url_file.write("\n")
        print("Microreact: " + r.url + "\n")
elif r.status_code == 400:
    sys.stderr.write("Microreact API call failed with error " + r.error + "\n")
else:
    sys.stderr.write("Microreact API call failed with unknown response code " + str(r.status_code) + "\n")
    sys.exit(1)
    