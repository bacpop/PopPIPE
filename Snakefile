import pandas as pd
from collections import defaultdict
from snakemake.utils import validate, min_version
min_version("5.3.0")

configfile: "config.yml"
validate(config, schema="config.schema.yml")

# table with sample, sequence, cluster
samples = pd.read_table(config["poppunk_rfile"], header=None, index_col=0)
clusters = pd.read_table(config["poppunk_clusters"], sep=",",).set_index("Taxon")

included_strain_ids = list((clusters.Cluster.value_counts()[clusters.Cluster.value_counts() >= config["min_cluster_size"]]).index)
included_samples = samples.loc[clusters.index[clusters.isin(included_strain_ids)["Cluster"]]]

# First rule is the default target
rule cluster_summary:
    input:
        expand("output/strains/{strain}/fastbaps_clusters.txt", strain=included_strain_ids)
    output:
        "output/all_clusters.txt"
    params:
        levels=config['fastbaps']['levels']
    script:
       "{config[poppipe_location]}/scripts/number_clusters.py"

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

rule group_stragglers:
    input:
        config["poppunk_rfile"]
    output:
        rfile="output/strains/other/rfile.txt",
        namefile="output/strains/other/names.txt"
    group:
        "clustersplit"
    run:
        excluded_strains = clusters.loc[clusters.index[~clusters.isin(included_strain_ids)["Cluster"]]]
        d = defaultdict(list)
        for strain in included_strain_ids:
            strain_list = clusters.loc[clusters['Cluster'] == int(strain)]
            d['Taxon'].append(strain_list.index[0])
            d['Cluster'].append(strain_list['Cluster'][0])
        all_df = pd.concat([pd.DataFrame(data=d).set_index('Taxon'), excluded_strains])
        all_df.to_csv(output.rfile, sep="\t", header=False)
        all_df.index.to_series().to_csv(output.namefile, sep="\t", header=False, index=False)

# Use sketchlib to extract distances
# TODO write wrapper script which generates(/updates) input from a PopPUNK run
# needs to delete output files from strain directories
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
        config["poppipe_location"] + "/envs/sketch.yml"
    shell:
        # TODO use this when I have fixed the subset option on the conda-forge version
        #"poppunk_sketch --query --ref-db {params.db_prefix} --query-db {params.db_prefix} --subset {input.names} "
        #"--output {params.dist_prefix} --min-k {params.min_k} --max-k {params.max_k} --k-step {params.k_step} "
        #"--sketch-size {params.sketch_size} --cpus {threads} &> {log}"
        "export DYLD_FALLBACK_LIBRARY_PATH=/Users/jlees/miniconda3/envs/gfortran/lib && python ../pp-sketchlib/pp_sketch-runner.py --query --ref-db {params.db_prefix} --query-db {params.db_prefix} "
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
        config["poppipe_location"] + "/envs/nj.yml"
    script:
        "{config[poppipe_location]}/scripts/run_rapidnj.py"

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
        config["poppipe_location"] + "/envs/ska.yml"
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
        config["poppipe_location"] + "/envs/ska.yml"
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
        config["poppipe_location"] + "/envs/iqtree.yml"
    threads:
        4
    script:
        "{config[poppipe_location]}/scripts/run_iqtree.py"

#TODO replace with Nick's lineage clustering algorithm
rule lineage_clust:
    input:
        npy="output/strains/{strain}/dists.npy",
        pkl="output/strains/{strain}/dists.pkl"
    output:
        "output/strains/{strain}/hclust.txt"
    group:
        "quickclust"
    script:
       "{config[poppipe_location]}/scripts/run_hclust.py" #TODO call the cluster numbering script from within here

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
    conda:
        config["poppipe_location"] + "/envs/fastbaps.yml"
    threads:
        2
    shell:
        "{params.fb_script} -p 'baps' -i {input.align} -o {output} -l {params.levels} --phylogney={input.tree} -t {threads} > {log}"

# Visualisation below here:

rule graft_tree:
    input:
        ml_trees=expand("output/strains/{strain}/besttree.nwk", strain=included_strain_ids),
        nj_trees=expand("output/strains/{strain}/njtree.nwk", strain=included_strain_ids),
        overall_tree="output/strains/other/njtree.nwk"
    output:
        "output/full_tree.nwk"
    group:
        "viz"
    conda:
        config["poppipe_location"] + "/envs/nj.yml"
    script:
        "{config[poppipe_location]}/scripts/tree_graft.py"

rule generate_dot:
    input:
        database = config["poppunk_db"] + ".h5"
    output:
        "output/all_dists.npy",
        "output/all_dists.pkl",
        tsne_out="output/embedding.dot"
    group:
        "viz"
    params:
        dist_prefix = "output/all_dists",
        db_prefix = config["poppunk_db"],
        perplexity = config['tsne']['perplexity']
    threads:
        16
    log:
        sketch_log = "logs/all_query.log",
        tsne_log = "logs/tsne.log"
    conda:
        config["poppipe_location"] + "/envs/poppunk.yml"
    shell:
        """
        export DYLD_FALLBACK_LIBRARY_PATH=/Users/jlees/miniconda3/envs/gfortran/lib && python ../pp-sketchlib/pp_sketch-runner.py --query --ref-db {params.db_prefix} --query-db {params.db_prefix} --output {params.dist_prefix} --read-k --cpus {threads} &> {log.sketch_log}
        poppunk_tsne --distances {params.dist_prefix} --output output --perplexity {params.perplexity} --verbosity 1 &> {log.tsne_log}
        mv output/output_perplexity{params.perplexity}_accessory_tsne.dot {output.tsne_out}
        """
        #"""
        #poppunk_sketch --query --ref-db {params.db_prefix} --query-db {params.db_prefix} --output {params.dist_prefix} --read-k --cpus {threads} &> {log.sketch_log}
        #poppunk_tsne --distances {params.dist_prefix} --output output --perplexity {params.perplexity} --verbosity 1 &> {log.tsne_log}
        #mv output/output_perplexity{params.perplexity}_accessory_tsne.dot {output.tsne_out}
        #"""

# use microreact api
rule make_microreact:
    input:
        tree="output/full_tree.nwk",
        clusters="output/all_clusters.txt",
        dot="output/embedding.dot"
    output:
        "output/microreact_url.txt"
    params:
        microreact_name=config['microreact']['name'],
        microreact_email=config['microreact']['email'],
        microreact_website=config['microreact']['website']
    group:
        "viz"
    script:
        "{config[poppipe_location]}/scripts/make_microreact.py"

#TODO possible extension - use snap for when the input is fastq
# snap for alignment - pick a good N50 to align to
rule snap_index:

# need to map, call and convert vcf to msa
rule snap_map:






