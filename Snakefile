import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.1.2")

configfile: "config.yml"
validate(config, schema="config.schema.yml")

# table with sample, sequence, cluster
samples = pd.read_table(config["poppunk_rfile"], header=None, index_col=0)
clusters = pd.read_table(config["poppunk_clusters"], sep=",",).set_index("Taxon")

pruned_clusters = list((clusters.Cluster.value_counts()[clusters.Cluster.value_counts() >= config["min_cluster_size"]]).index)
pruned_samples = samples.loc[clusters.index[clusters.isin(pruned_clusters)["Cluster"]]]

#TODO python script to generate config.yml (use yaml module)
#this could also execute snakemake targeting specific rules, if certain steps are being skipped for speed

# First rule is the default target
rule all:
    input:
        "all_clusters.txt"

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
        "logs/{cluster}_sketchlib.log"
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
        assembly=lambda wildcards: pruned_samples.loc[wildcards.sample]
    output:
        "ska_index/{sample}.skf"
    log:
        "logs/ska_index_{sample}.log" 
    conda:
        "envs/ska.yml"
    shell:
        "ska fasta -o ska_index/{wildcards.sample} {input.assembly} > {log}"

rule ska_align:
    input:
        ska=expand("ska_index/{sample}.skf", sample=pruned_samples.index),
        samples="clusters/{cluster}/ska_index.txt"
    output:
        "clusters/{cluster}/align_variants.aln"
    log:
        "logs/ska_align_{cluster}.log" 
    params:
        prefix="clusters/{cluster}/align"
    conda:
        "envs/ska.yml"
    shell:
        "ska align -v -o {params.prefix} -f {input.samples} > {log}"

# will set fast or slow in params.yaml
rule iq_tree:
    input:
        start_tree="clusters/{cluster}/njtree.nwk",
        alignment="clusters/{cluster}/align_variants.aln"
    output:
        rooted="clusters/{cluster}/besttree.nwk",
        unrooted=temp("clusters/{cluster}/besttree.unrooted.treefile"),
        iqtree="clusters/{cluster}/besttree.unrooted.iqtree",
        ckp="clusters/{cluster}/besttree.unrooted.ckp.gz",
    log:
        "clusters/{cluster}/besttree.unrooted.log"
    params:
        enabled=config['iqtree']['enabled'],
        mode=config['iqtree']['mode'],
        model=config['iqtree']['model'],
        prefix="clusters/{cluster}/besttree.unrooted"
    conda:
        "envs/iqtree.yml"
    threads:
        4
    script:
        "{config[script_location]}/run_iqtree.py"

# in tree + aln mode
rule fastbaps:
    input:
        tree="clusters/{cluster}/besttree.nwk",
        align="clusters/{cluster}/align_variants.aln"
    output:
        "clusters/{cluster}/fastbaps_clusters.txt"
    params:
        fb_script=config['fastbaps']['script'],
        levels=config['fastbaps']['levels']
    log:
        "logs/fastbaps_{cluster}.log"
    #conda:
    #    "envs/fastbaps.yml"
    threads:
        4
    shell:
        "{params.fb_script} -i {input.align} -o {output} -l {params.levels} --phylogney={input.tree} -t {threads} > {log}"
        
rule cluster_summary:
    input:
       expand("clusters/{cluster}/fastbaps_clusters.txt", cluster=pruned_clusters),
    output:
       "all_clusters.txt"
    shell:
       "cat {input} > {output}"


# snap for alignment
rule snap_index:

rule snap_map:
# need to convert vcf to msa



# run overall rapidnj and t-sne
rule generate_viz:

# use microreact api
rule make_microreact:


