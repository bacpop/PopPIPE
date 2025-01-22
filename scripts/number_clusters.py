import os, sys
import re
import pandas as pd

samples = pd.read_table(snakemake.config["poppunk_rfile"], header=None, index_col=0)
clusters = pd.read_table(snakemake.config["poppunk_clusters"], sep=",",).set_index("Taxon")

included_strain_ids = list((clusters.Cluster.value_counts()[clusters.Cluster.value_counts() >= snakemake.config["min_cluster_size"]]).index)
excluded_clusters = clusters.loc[clusters.index[~clusters.isin(included_strain_ids)["Cluster"]]]

with open(snakemake.output[0], 'w') as outfile:
    outfile.write(",".join(["Taxon", "Strain"] + \
                           ["Subcluster_" + str(x + 1) for x in range(snakemake.params['levels'])])
                  + "\n")
    cluster_max = [0] * snakemake.params['levels']

    # Read each fastbaps file
    for fb_file in snakemake.input:
        # Get PopPUNK cluster/strain number
        file_match = re.search(r'output\/strains\/([\w_]+)\/fastbaps_clusters.txt$', fb_file)
        if file_match:
            cluster = file_match.group(1)
        else:
            print("no match: ", file_match)
            sys.stderr.write("Error finding fastbaps output\n")
            sys.exit(1)

        # Read the fastbaps clusters
        with open(fb_file, 'r') as fastbaps_clusters:
            header = fastbaps_clusters.readline()
            prev_max = cluster_max.copy()

            for cluster_line in fastbaps_clusters:
                fields_in = cluster_line.rstrip().split(",")
                fields_out = [fields_in[0], cluster]

                # For each level, number clusters consecutively
                for subcluster_idx, subcluster in enumerate(fields_in[1:]):
                    new_subcluster = int(subcluster) + prev_max[subcluster_idx]
                    fields_out += [str(new_subcluster)]
                    if new_subcluster > cluster_max[subcluster_idx]:
                        cluster_max[subcluster_idx] = new_subcluster
                outfile.write(",".join(fields_out) + "\n")

    for sample, cluster in excluded_clusters.iterrows():
        fields_out = [sample, str(cluster['Cluster'])] + ['NA'] * snakemake.params['levels']
        outfile.write(",".join(fields_out) + "\n")
