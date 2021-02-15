#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :

import sys, os
import pandas as pd
import re
import subprocess

import yaml

# script outputs
combined_rfile = "combined_rfile.txt"
combined_clusters = "combined_clusters.csv"
combined_db = "combined_db.h5"

# command line parsing
def get_options():
    import argparse
    parser = argparse.ArgumentParser(description='Update PopPIPE results with queries',
                                     prog='poppipe_assign.py')

    # input options
    qGroup = parser.add_argument_group('poppunk_assign files')
    qGroup.add_argument('--query', required=True, help='Argument given to --query')
    qGroup.add_argument('--db', required=True, help='Argument given to --db')
    qGroup.add_argument('--output', required=True, help='Argument given to --output')

    smGroup = parser.add_argument_group('snakemake files')
    smGroup.add_argument('--rfile', required=True, help='poppunk_rfile from previous run')
    smGroup.add_argument('--clusters', required=True, help='poppunk_clusters from previous run')
    smGroup.add_argument('--h5', required=True, help='poppunk_h5 from previous run')

    other = parser.add_argument_group('Other options')
    smGroup.add_argument('--recalculate-random', default=False, action='store_true',
                          help='Recalculate random matches when joining databases')

    return parser.parse_args()

# main code
if __name__ == "__main__":

    # Get command line options
    args = get_options()

    # Join files
    with open(combined_rfile, 'w') as rfile_out:
        with open(args.rfile, 'r') as r1_in:
            rfile_out.write(r1_in.read())
        with open(args.query, 'r') as r2_in:
            rfile_out.write(r2_in.read())

    # Join clusters
    ref_clusters = pd.read_table(args.clusters, sep=",",).set_index("Taxon")
    query_clusters = pd.read_table(
        args.output + "/" + os.path.basename(args.output) + "_clusters.csv",
        sep=",",).set_index("Taxon")
    ref_clusters.append(query_clusters)
    pd.to_csv(combined_clusters, sep=",")

    # Join DBs
    h5_re = re.compile(r"^(.+)\.h5$")
    ref_prefix = re.match(h5_re, args.h5).group(1)
    output_prefix = re.match(h5_re, ref_prefix).group(1)
    query_prefix = args.output + "/" + os.path.basename(args.output)

    join_cmd = "poppunk_sketch --join --ref-db " + ref_prefix + \
               " --query-db " + query_prefix + " --output " + \
               output_prefix
    subprocess.run(join_cmd, shell=True, check=True)
    if args.recalculate_random:
        subprocess.run("poppunk_sketch --add-random --ref-db" + \
                        output_prefix,
                       shell=True, check=True)

    with open("config.yaml", 'r') as stream:
        try:
            config_in = yaml.safe_load(stream))
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit(1)

    config_in['poppunk_rfile'] = combined_rfile
    config_in['poppunk_clusters'] = combined_clusters
    config_in['poppunk_h5'] = combined_db

    with open("config.yaml.test", 'w') as stream:
        output = yaml.dump(config_in)

    sys.exit(0)
