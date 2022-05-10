import sys
import pandas as pd
from datetime import datetime
import pickle
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
description_string = "PopPIPE run on " + datetime.now().strftime("%Y-%b-%d %H:%M")

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

# Load example JSON to be modified
with open(snakemake.params['example_file'], 'rb') as example_pickle:
    json_pickle = pickle.load(example_pickle)

json_pickle["files"]["data-file-1"]["blob"] = csv_string
json_pickle["files"]["tree-file-1"]["blob"] = tree_string
json_pickle["files"]["network-file-1"] = {"id": "network-file-1",
                                          "name": "network.dot",
                                          "format": "text/vnd.graphviz",
                                          "blob": dot_string}
json_pickle["networks"]["network-1"] = {"title": "Network",
                                        "file": "network-file-1",
                                        "nodeField": "id"}
json_pickle["meta"]["name"] = snakemake.params['microreact_name']
json_pickle["meta"]["description"] = description_string
json_pickle["meta"]["email"] = snakemake.params['microreact_email']
json_pickle["meta"]["email"] = snakemake.params['microreact_website']


# Use API to convert to new JSON
headers = {"Content-type": "application/json; charset=UTF-8",
           "Access-Token": api_token}
create_request = make_request(json.dumps(json_pickle), microreact_api_new_url, headers)

with open(snakemake.output['microreact'], 'w') as json_file:
    json_file.write(create_request.text)
    json_file.write("\n")

with open(snakemake.output['url'], 'w') as url_file:
    url = create_request.json()['url']
    url_file.write(url)
    url_file.write("\n")
    print("Microreact: " + url + "\n")

