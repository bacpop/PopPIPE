import os, sys
import re

with open(snakemake.output[0], 'w') as outfile:
    outfile.write(",".join(["Taxon", "Strain"] + \
                           ["Subcluster_" + str(x + 1) for x in range(snakemake.params['levels'])])
                  + "\n")
    cluster_max = [0] * snakemake.params['levels'] 
    
    # Read each fastbaps file
    for fb_file in snakemake.input:
        # Get PopPUNK cluster/strain number
        file_match = re.search(r'^output\/strains\/(\d+)\/fastbaps_clusters.txt$', fb_file)
        if file_match:
            cluster = file_match.group(1)
        else:
            sys.stderr.write("Error finding fastbaps output\n")
            sys.exit(1)
        with open(fb_file, 'r') as fastbaps_clusters:
            header = fastbaps_clusters.readline()
            prev_max = cluster_max.copy()
            
            for cluster_line in fastbaps_clusters:
                fields_in = cluster_line.rstrip().split(",")
                fields_out = [fields_in[0], cluster]
                for subcluster_idx, subcluster in enumerate(fields_in[1:]):
                    new_subcluster = int(subcluster) + prev_max[subcluster_idx]
                    fields_out += [str(new_subcluster)]
                    if new_subcluster > cluster_max[subcluster_idx]:
                        cluster_max[subcluster_idx] = new_subcluster
                outfile.write(",".join(fields_out) + "\n")
                    


