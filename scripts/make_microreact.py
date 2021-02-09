import os, sys
import pandas as pd
from datetime import date
import requests

# Constants
microreact_api_url = "https://microreact.org/api/project/"
description_string = "PopPIPE run on " + date.today().strftime("%Y-%b-%d %H:%M")

# Format data
with open(snakemake.input['tree'], 'r') as tree_file:
    tree_string = tree_file.read()

with open(snakemake.input['dot'], 'r') as dot_file:
    dot_string = dot_file.read()

clusters = pd.read_table(snakemake.input['clusters'], sep=",", dtype=str,
                         keep_default_na=False, na_values=[""]).set_index("Taxon")
clusters.index.rename('id', inplace=True)
clusters.rename(columns={i: i + "__autocolour" for i in clusters.columns}, inplace=True)
csv_string = clusters.to_csv(None, sep=",", header=True, index=True)

# Make request
payload = {"name": snakemake.params['microreact_name'],
           "description": description_string,
           "email": snakemake.params['microreact_email'],
           "website": snakemake.params['microreact_website'],
           "data": csv_string,
           "tree": tree_string,
           "dot": dot_string}

r = requests.post(microreact_api_url, data = payload)
if r.ok:
    url = r.json()['url']
    with open(snakemake.output[0], 'w') as url_file:
        url_file.write(url)
        url_file.write("\n")
        print("Microreact: " + url + "\n")
elif r.status_code == 400:
    sys.stderr.write("Microreact API call failed with response " + r.json()['error'] + "\n")
else:
    sys.stderr.write("Microreact API call failed with unknown response code " + str(r.status_code) + "\n")
    sys.exit(1)
