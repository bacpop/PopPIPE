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

# First rule is the default target
rule all:
    input:
        "output/all_clusters.txt"

rule split_clusters:
    input:
        config["poppunk_rfile"]
    output:
        rfile="output/clusters/{cluster}/rfile.txt",
        namefile="output/clusters/{cluster}/names.txt",
        skafile="output/clusters/{cluster}/ska_index.txt"
    group:
        "clustersplit"
    run:
        cluster_samples = clusters.loc[clusters['Cluster'] == int(wildcards.cluster)].index
        sample_subset = samples.loc[list(cluster_samples)]
        if (len(sample_subset) > 1):
            sample_subset.to_csv(output.rfile, sep="\t", header=False)
            sample_subset.index.to_series().to_csv(output.namefile, sep="\t", header=False, index=False)
            with open(output.skafile, 'w') as ska_list:
                for sample in cluster_samples:
                    ska_list.write("output/ska_index/" + sample + ".skf\n")

# Use sketchlib to extract distances
# TODO write wrapper script which generates(/updates) input from a PopPUNK run
rule sketchlib_dists:
    input:
        database = config["poppunk_db"] + ".h5",
        names = "output/clusters/{cluster}/names.txt"
    output:
        "output/clusters/{cluster}/dists.npy",
        "output/clusters/{cluster}/dists.pkl"
    group:
        "sketchlib"
    log:
        "logs/{cluster}_sketchlib.log"
    params:
        db_prefix = config["poppunk_db"],
        dist_prefix = "output/clusters/{cluster}/dists"
    threads:
        16
    conda:
        "envs/sketch.yml"
    shell:
        # TODO use this when I have fixed the subset option on the conda-forge version
        #"poppunk_sketch --query --ref-db {params.db_prefix} --query-db {params.db_prefix} --subset {input.names} "
        #"--output {params.dist_prefix} --min-k {params.min_k} --max-k {params.max_k} --k-step {params.k_step} "
        #"--sketch-size {params.sketch_size} --cpus {threads} &> {log}"
        "python ../pp-sketchlib/pp_sketch-runner.py --query --ref-db {params.db_prefix} --query-db {params.db_prefix} "
        "--subset {input.names} --output {params.dist_prefix} --read-k --cpus {threads} &> {log}"

# rapidnj
rule generate_nj:
    input:
        npy="output/clusters/{cluster}/dists.npy",
        pkl="output/clusters/{cluster}/dists.pkl" 
    output:
        "output/clusters/{cluster}/njtree.nwk"
    group:
        "quicktree"
    conda:
        "envs/nj.yml"
    script:
        "{config[script_location]}/run_rapidnj.py"

# ska for alignment
#TODO modify this to work out whether fastq version is needed
rule ska_index:
    input:
        assembly=lambda wildcards: pruned_samples.loc[wildcards.sample]
    output:
        "output/ska_index/{sample}.skf"
    log:
        "logs/ska_index_{sample}.log" 
    conda:
        "envs/ska.yml"
    shell:
        "ska fasta -o output/ska_index/{wildcards.sample} {input.assembly} > {log}"

rule ska_align:
    input:
        ska=expand("output/ska_index/{sample}.skf", sample=pruned_samples.index),
        samples="output/clusters/{cluster}/ska_index.txt"
    output:
        "output/clusters/{cluster}/align_variants.aln"
    group:
        "align"
    log:
        "logs/ska_align_{cluster}.log" 
    params:
        prefix="output/clusters/{cluster}/align"
    conda:
        "envs/ska.yml"
    shell:
        "ska align -v -o {params.prefix} -f {input.samples} > {log}"

# will set fast or slow in params.yaml
rule iq_tree:
    input:
        start_tree="output/clusters/{cluster}/njtree.nwk",
        alignment="output/clusters/{cluster}/align_variants.aln"
    output:
        rooted="output/clusters/{cluster}/besttree.nwk",
        unrooted=temp("output/clusters/{cluster}/besttree.unrooted.treefile"),
        iqtree=temp("output/clusters/{cluster}/besttree.unrooted.iqtree"),
        ckp=temp("output/clusters/{cluster}/besttree.unrooted.ckp.gz")
    log:
        "output/clusters/{cluster}/besttree.unrooted.log"
    params:
        enabled=config['iqtree']['enabled'],
        mode=config['iqtree']['mode'],
        model=config['iqtree']['model'],
        prefix="output/clusters/{cluster}/besttree.unrooted"
    conda:
        "envs/iqtree.yml"
    threads:
        4
    script:
        "{config[script_location]}/run_iqtree.py"

rule hclust:
    input:
        npy="output/clusters/{cluster}/dists.npy",
        pkl="output/clusters/{cluster}/dists.pkl" 
    output:
        "output/clusters/{cluster}/hclust.txt"
    group:
        "quickclust"
    script:
       "{config[script_location]}/run_hclust.py" #TODO call the cluster numbering script from within here

# in tree + aln mode
rule fastbaps:
    input:
        tree="output/clusters/{cluster}/besttree.nwk",
        align="output/clusters/{cluster}/align_variants.aln"
    output:
        "output/clusters/{cluster}/fastbaps_clusters.txt"
    group:
        "bapsclust"
    params:
        fb_script=config['fastbaps']['script'],
        levels=config['fastbaps']['levels']
    log:
        "logs/fastbaps_{cluster}.log"
    #conda:
    #    "envs/fastbaps.yml"
    threads:
        2
    shell:
        "{params.fb_script} -p 'baps' -i {input.align} -o {output} -l {params.levels} --phylogney={input.tree} -t {threads} > {log}"
        
rule cluster_summary:
    input:
       expand("output/clusters/{cluster}/fastbaps_clusters.txt", cluster=pruned_clusters),
    output:
       "output/all_clusters.txt"
    shell:
       "cat {input} > {output}" #TODO replace this with a script to give secondary clusters consistent names


# run overall rapidnj and t-sne
rule generate_viz:

# use microreact api
rule make_microreact:

#TODO possible extension - use snap for when the input is fastq
# snap for alignment - pick a good N50 to align to
rule snap_index:

# need to map, call and convert vcf to msa
rule snap_map:






