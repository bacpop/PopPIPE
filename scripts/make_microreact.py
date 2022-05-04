import os, sys
import pandas as pd
from datetime import date
import requests
import json

def make_request(payload, url, headers):
    r = requests.post(url, data=payload, headers=headers)

    if not r.ok:
        if r.status_code == 400:
            sys.stderr.write("Microreact API call failed with response " + r.text + "\n")
        else:
            sys.stderr.write("Microreact API call failed with unknown response code " + str(r.status_code) + "\n")
        sys.exit(1)
    return r

# Constants
api_token = snakemake.params['microreact_token']
microreact_api_convert_url = "https://microreact.org/api/schema/convert"
microreact_api_new_url = "https://microreact.org/api/projects/create"
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

# Use API to convert to new JSON
headers = {"Content-type": "application/json; charset=UTF-8"}
payload = {"name": snakemake.params['microreact_name'],
           "description": description_string,
           "email": snakemake.params['microreact_email'],
           "website": snakemake.params['microreact_website'],
           "data": csv_string,
           "tree": tree_string,
           "dot": dot_string}
new_json = make_request(json.dumps(payload), microreact_api_convert_url, headers)

# Use this to create new microreact
headers['Access-Token'] = api_token
create_request = make_request(new_json.text, microreact_api_new_url, headers)

url = create_request.json()['url']
with open(snakemake.output[0], 'w') as url_file:
    url_file.write(url)
    url_file.write("\n")
    print("Microreact: " + url + "\n")

