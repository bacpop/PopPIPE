import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.1.2")

configfile: "config.yml"
validate(config, schema="config.schema.yml")

# table with sample, sequence, cluster
samples = pd.read_table(config["poppunk_rfile"], header=None, index_col=0)
clusters = pd.read_table(config["poppunk_clusters"], sep=",",).set_index("Taxon", drop=False)
cluster_names=pd.unique(clusters['Cluster'])
# all or a subset?
# ideally automatically update an existing analysis with new isolates
# general, python script to generate a params.yaml? 
# Or just edit a template

# First rule is the default target
rule all:
    input:
        "output/cluster_summary.txt"

rule split_clusters:
    input:
        config["poppunk_rfile"]
    output:
        "clusters/{cluster}/rfile.txt",
        "clusters/{cluster}/ska_index.txt"
    script:
        "{config[script_loc]}/update_cluster_dirs.py" #TODO script to read PopPUNK results, create cluster dirs, add any new samples

# Use sketchlib to extract distances
# TODO write wrapper script which generates(/updates) input from a PopPUNK run
rule sketchlib_dists:
    input:
        database = config["poppunk_db"],
        rfile = "clusters/{cluster}/rfile.txt"
    output:
        "clusters/{cluster}/dists.npy",
        "clusters/{cluster}/dists.pkl"
    conda:
        "envs/sketch.yml"
    shell:
        "poppunk_sketch --query --ref-db {input.database} --query-db {input.database} --subset {input.rfile} "
        "--output clusters/{cluster}/dists --min-k {config[min_k]} --max-k {config[max_k]} --k-step {config[k_step]} "
        "--sketch-size {config[sketch_size]} --cpus {threads}"

# rapidnj
rule generate_nj:
    input:
        npy="clusters/{cluster}/dists.npy",
        pkl="clusters/{cluster}/dists.pkl" 
    output:
        "clusters/{cluster}/njtree.nwk"
    conda:
        "envs/nj.yml"
    script:
        "{config[script_location]}/run_rapidnj.py"  #TODO use function from PopPUNK

# ska for alignment
rule ska_index:
    input:
        assembly=lambda wildcards: samples.loc[wildcards.sample]
    output:
        "ska_index/{sample}.skf"
    conda:
        "envs/ska.yml"
    shell:
        "ska fasta -o ska_index/{wildcards.sample} {input.assembly}"

rule ska_align:
    input:
        ska=expand("ska_index/{sample}.skf", sample=clusters.index),
        samples="clusters/{cluster}/ska_index.txt"
    output:
        "clusters/{cluster}/align.fa",
    conda:
        "envs/ska.yml"
    shell:
        "ska align -v -o clusters/{cluster}/align -f {input.samples}"

#TODO remove this
rule dummy_finish:
    input:
       expand("clusters/{cluster}/align.fa", cluster=cluster_names),
       expand("clusters/{cluster}/njtree.nwk", cluster=cluster_names)
    output:
       "output/cluster_summary.txt"
    shell:
        "touch {output}"
 


# first a rule which cats ska file names into new list
# then run align using this file

### GET TO HERE WORKING BEFORE PROCEEDING BELOW

# snap for alignment
rule snap_index:

rule snap_map:
# need to convert vcf to msa

# will set fast or slow in params.yaml
rule iq_tree:

# in tree + aln mode
rule fastbaps:

# run overall rapidnj and t-sne
rule generate_viz:

# use microreact api
rule make_microreact:


