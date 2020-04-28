import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.1.2")

configfile: "config.yml"
validate(config, schema="config.schema.yml")

# table with sample, sequence, cluster
samples = pd.read_table(config["poppunk_rfile"], header=None, index_col=0)
clusters = pd.read_table(config["poppunk_clusters"], sep=",",).set_index("Taxon")

included_strain_ids = list((clusters.Cluster.value_counts()[clusters.Cluster.value_counts() >= config["min_cluster_size"]]).index)
included_samples = samples.loc[clusters.index[clusters.isin(included_strain_ids)["Cluster"]]]

# First rule is the default target
rule all:
    input:
        "output/all_clusters.txt"

rule split_strains:
    input:
        config["poppunk_rfile"]
    output:
        rfile="output/strains/{strain}/rfile.txt",
        namefile="output/strains/{strain}/names.txt",
        skafile="output/strains/{strain}/ska_index.txt"
    group:
        "clustersplit"
    run:
        cluster_samples = clusters.loc[clusters['Cluster'] == int(wildcards.strain)].index
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
        names = "output/strains/{strain}/names.txt"
    output:
        "output/strains/{strain}/dists.npy",
        "output/strains/{strain}/dists.pkl"
    group:
        "sketchlib"
    log:
        "logs/{strain}_sketchlib.log"
    params:
        db_prefix = config["poppunk_db"],
        dist_prefix = "output/strains/{strain}/dists"
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
        npy="output/strains/{strain}/dists.npy",
        pkl="output/strains/{strain}/dists.pkl" 
    output:
        "output/strains/{strain}/njtree.nwk"
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
        assembly=lambda wildcards: included_samples.loc[wildcards.sample]
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
        ska=expand("output/ska_index/{sample}.skf", sample=included_samples.index),
        samples="output/strains/{strain}/ska_index.txt"
    output:
        "output/strains/{strain}/align_variants.aln"
    group:
        "align"
    log:
        "logs/ska_align_{strain}.log" 
    params:
        prefix="output/strains/{strain}/align"
    conda:
        "envs/ska.yml"
    shell:
        "ska align -v -o {params.prefix} -f {input.samples} > {log}"

# will set fast or slow in params.yaml
rule iq_tree:
    input:
        start_tree="output/strains/{strain}/njtree.nwk",
        alignment="output/strains/{strain}/align_variants.aln"
    output:
        rooted="output/strains/{strain}/besttree.nwk",
        unrooted=temp("output/strains/{strain}/besttree.unrooted.treefile"),
        iqtree=temp("output/strains/{strain}/besttree.unrooted.iqtree"),
        ckp=temp("output/strains/{strain}/besttree.unrooted.ckp.gz")
    log:
        "output/strains/{strain}/besttree.unrooted.log"
    params:
        enabled=config['iqtree']['enabled'],
        mode=config['iqtree']['mode'],
        model=config['iqtree']['model'],
        prefix="output/strains/{strain}/besttree.unrooted"
    conda:
        "envs/iqtree.yml"
    threads:
        4
    script:
        "{config[script_location]}/run_iqtree.py"

rule hclust:
    input:
        npy="output/strains/{strain}/dists.npy",
        pkl="output/strains/{strain}/dists.pkl" 
    output:
        "output/strains/{strain}/hclust.txt"
    group:
        "quickclust"
    script:
       "{config[script_location]}/run_hclust.py" #TODO call the cluster numbering script from within here

# in tree + aln mode
rule fastbaps:
    input:
        tree="output/strains/{strain}/besttree.nwk",
        align="output/strains/{strain}/align_variants.aln"
    output:
        "output/strains/{strain}/fastbaps_clusters.txt"
    group:
        "bapsclust"
    params:
        fb_script=config['fastbaps']['script'],
        levels=config['fastbaps']['levels']
    log:
        "logs/fastbaps_{strain}.log"
    #conda:
    #    "envs/fastbaps.yml"
    threads:
        2
    shell:
        "{params.fb_script} -p 'baps' -i {input.align} -o {output} -l {params.levels} --phylogney={input.tree} -t {threads} > {log}"
        
rule cluster_summary:
    input:
        expand("output/strains/{strain}/fastbaps_clusters.txt", strain=included_strain_ids)
    output:
        "output/all_clusters.txt"
    params:
        levels=config['fastbaps']['levels']
    script:
       "{config[script_location]}/number_clusters.py"


# run overall rapidnj and t-sne
rule generate_viz:
#TODO refine overall tree by replacing subclades with ML trees

# use microreact api
rule make_microreact:

#TODO possible extension - use snap for when the input is fastq
# snap for alignment - pick a good N50 to align to
rule snap_index:

# need to map, call and convert vcf to msa
rule snap_map:






