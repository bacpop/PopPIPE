import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.1.2")

configfile: "config.yml"
validate(config, schema="config.schema.yml")

# table with sample, sequence, cluster
samples = pd.read_table(config["poppunk_rfile"], header=None, index_col=0)
clusters = pd.read_table(config["poppunk_clusters"], sep=",",).set_index("Taxon")
non_singleton_clusters = list((clusters.Cluster.value_counts()[clusters.Cluster.value_counts() >= config["min_cluster_size"]]).index)

#TODO only index ska files used in clusters

#TODO use dynamic() to make ska files specific to cluster

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
        rfile="clusters/{cluster}/rfile.txt",
        namefile="clusters/{cluster}/names.txt",
        skafile="clusters/{cluster}/ska_index.txt"
    run:
        cluster_samples = clusters.loc[clusters['Cluster'] == int(wildcards.cluster)].index
        sample_subset = samples.loc[list(cluster_samples)]
        if (len(sample_subset) > 1):
            sample_subset.to_csv(output.rfile, sep="\t", header=False)
            sample_subset.index.to_series().to_csv(output.namefile, sep="\t", header=False, index=False)
            with open(output.skafile, 'w') as ska_list:
                for sample in cluster_samples:
                    ska_list.write("ska_index/" + sample + ".skf\n")

# Use sketchlib to extract distances
# TODO write wrapper script which generates(/updates) input from a PopPUNK run
rule sketchlib_dists:
    input:
        database = config["poppunk_db"] + ".h5",
        names = "clusters/{cluster}/names.txt"
    output:
        "clusters/{cluster}/dists.npy",
        "clusters/{cluster}/dists.pkl"
    log:
        "clusters/{cluster}/sketchlib.log"
    params:
        db_prefix = config["poppunk_db"],
        dist_prefix = "clusters/{cluster}/dists",
        min_k = config['sketch']['min_k'],
        max_k = config['sketch']['max_k'],
        k_step = config['sketch']['k_step'],
        sketch_size = config['sketch']['sketch_size']
    threads:
        16
    conda:
        "envs/sketch.yml"
    shell:
        # TODO use this when I have fixed the subset option on the conda-forge version
        #"poppunk_sketch --query --ref-db {params.db_prefix} --query-db {params.db_prefix} --subset {input.names} "
        #"--output {params.dist_prefix} --min-k {params.min_k} --max-k {params.max_k} --k-step {params.k_step} "
        #"--sketch-size {params.sketch_size} --cpus {threads} &> {log}"
        "python ../pp-sketchlib/pp_sketch-runner.py --query --ref-db {params.db_prefix} --query-db {params.db_prefix} --subset {input.names} "
        "--output {params.dist_prefix} --min-k {params.min_k} --max-k {params.max_k} --k-step {params.k_step} "
        "--sketch-size {params.sketch_size} --cpus {threads} &> {log}"

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
        "{config[script_location]}/run_rapidnj.py"

# ska for alignment
rule ska_index:
    input:
        assembly=lambda wildcards: samples.loc[wildcards.sample]
    output:
        "ska_index/{sample}.skf"
    log:
        "ska_index/{sample}.log" 
    conda:
        "envs/ska.yml"
    shell:
        "ska fasta -o ska_index/{wildcards.sample} {input.assembly} > {log}"

rule ska_align:
    input:
        ska=expand("ska_index/{sample}.skf", sample=clusters.index),
        samples="clusters/{cluster}/ska_index.txt"
    output:
        "clusters/{cluster}/align_variants.aln"
    log:
       "clusters/{cluster}/ska.log" 
    params:
        prefix="clusters/{cluster}/align"
    conda:
        "envs/ska.yml"
    shell:
        "ska align -v -o {params.prefix} -f {input.samples} > {log}"

#TODO remove this
rule dummy_finish:
    input:
       expand("clusters/{cluster}/align_variants.aln", cluster=non_singleton_clusters),
       expand("clusters/{cluster}/njtree.nwk", cluster=non_singleton_clusters)
    output:
       "output/cluster_summary.txt"
    shell:
       "touch {output}"
 


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


