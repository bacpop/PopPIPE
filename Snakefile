
import re
import sys
from collections import defaultdict
import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.3.0")

configfile: "config.yml"
validate(config, schema="config.schema.yml")

# poppunk DB prefix
prefix_match = re.match(r"^(.+)\.h5$", config["poppunk_h5"])
if prefix_match:
    db_prefix = prefix_match.group(1)
else:
    raise RuntimeError("PopPUNK DB is not a .h5 file")

# table with sample, sequence, cluster
samples = pd.read_table(config["poppunk_rfile"], header=None, index_col=0)
clusters = pd.read_table(config["poppunk_clusters"], sep=",", dtype={'Cluster': str}).set_index("Taxon")

included_strain_ids = list((clusters.Cluster.value_counts()[clusters.Cluster.value_counts() >= config["min_cluster_size"]]).index)
if len(included_strain_ids) == 0:
    sys.exit("Error: No included strain IDs found. Check min_cluster_size.")
included_samples = samples.loc[clusters.index[clusters.isin(included_strain_ids)["Cluster"]]]

container: "docker://poppunk/poppipe:latest"

# First rule is the default target
rule cluster_summary:
    input:
        expand("output/strains/{strain}/fastbaps_clusters.txt", strain=included_strain_ids)
    output:
        "output/all_clusters.txt"
    params:
        levels=config['fastbaps']['levels']
    log:
        "logs/cluster_summary.log"
    script:
       config["poppipe_location"] + "/scripts/number_clusters.py"

rule split_strains:
    input:
        config["poppunk_rfile"]
    output:
        rfile="output/strains/{strain}/rfile.txt",
        namefile="output/strains/{strain}/names.txt"
    group:
        "clustersplit"
    run:
        cluster_samples = clusters.loc[clusters['Cluster'] == wildcards.strain].index
        sample_subset = samples.loc[list(cluster_samples)]
        if (len(sample_subset) > 1):
            sample_subset.to_csv(output.rfile, sep="\t", header=False)
            sample_subset.index.to_series().to_csv(output.namefile, sep="\t", header=False, index=False)

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
            strain = str(strain).strip()
            strain_list = clusters[clusters['Cluster'].astype(str).str.strip() == strain]
            if strain_list.empty:
                continue
            d['Taxon'].append(strain_list.index[0])
            d['Cluster'].append(strain_list['Cluster'].iloc[0])
        all_df = pd.concat([pd.DataFrame(data=d).set_index('Taxon'), excluded_strains])
        all_df.to_csv(output.rfile, sep="\t", header=False)
        all_df.index.to_series().to_csv(output.namefile, sep="\t", header=False, index=False)

# Use sketchlib to extract distances
rule sketchlib_dists:
    input:
        database = config["poppunk_h5"],
        names = "output/strains/{strain}/names.txt"
    output:
        "output/strains/{strain}/dists.npy",
        "output/strains/{strain}/dists.pkl"
    group:
        "sketchlib"
    log:
        "logs/{strain}_sketchlib.log"
    params:
        db_prefix = db_prefix,
        dist_prefix = "output/strains/{strain}/dists"
    threads:
        4
    shell:
        "sketchlib query dist {params.db_prefix} {params.db_prefix} --subset {input.names} "
        "-o {params.dist_prefix} --cpus {threads} &> {log}"

# rapidnj
rule generate_nj:
    input:
        npy="output/strains/{strain}/dists.npy",
        pkl="output/strains/{strain}/dists.pkl"
    output:
        start_tree="output/strains/{strain}/njtree.nwk"
    group:
        "quicktree"
    script:
        config["poppipe_location"] + "/scripts/run_rapidnj.py"

rule ska_build:
    input:
        samples="output/strains/{strain}/rfile.txt"
    output:
        skf="output/strains/{strain}/split_kmers.skf",
    params:
        skf_prefix="output/strains/{strain}/split_kmers",
        fastq_qual=config['ska']['fastq_cov'],
        fastq_cov=config['ska']['fastq_qual'],
        kmer=config['ska']['kmer'],
        single_strand=config['ska']['single_strand']
    log:
        "logs/ska_build_{strain}.log"
    script:
        config["poppipe_location"] + "/scripts/run_ska_build.py"

# ska for alignment
rule ska_align:
    input:
        skf="output/strains/{strain}/split_kmers.skf"
    output:
        alignment="output/strains/{strain}/align_variants.aln"
    group:
        "align"
    log:
        "logs/ska_align_{strain}.log"
    params:
        prefix="output/strains/{strain}/align"
    shell:
        "ska align -v --filter no-const --no-gap-only-sites {input.skf} > {output.alignment} 2> {log}"

# ska for mapping (needed for gubbins)
rule ska_map:
    input:
        skf="output/strains/{strain}/split_kmers.skf",
        reference_list="output/strains/{strain}/rfile.txt"
    output:
        alignment="output/strains/{strain}/map_variants.aln"
    group:
        "gubbins"
    log:
        "logs/ska_map_{strain}.log"
    params:
        prefix="output/strains/{strain}/align"
    shell:
        "ska map -v \"$(head -1 {input.reference_list} | cut -f 2)\" {input.skf} --ambig-mask > {output.alignment} 2> {log}"

# will set fast or slow in params.yaml
rule iq_tree:
    input:
        start_tree="output/strains/{strain}/njtree.nwk",
        alignment="output/strains/{strain}/align_variants.aln",
        rfiles=config["poppunk_rfile"]
    output:
        rooted="output/strains/{strain}/besttree.nwk",
        unrooted=temp("output/strains/{strain}/besttree.unrooted.treefile"),
        iqtree=temp("output/strains/{strain}/besttree.unrooted.iqtree")
    log:
        "output/strains/{strain}/besttree.unrooted.log"
    params:
        enabled=config['iqtree']['enabled'],
        mode=config['iqtree']['mode'],
        model=config['iqtree']['model'],
        prefix="output/strains/{strain}/besttree.unrooted"

    threads:
        4
    script:
        config["poppipe_location"] + "/scripts/run_iqtree.py"

rule gubbins:
    input:
        alignment="output/strains/{strain}/map_variants.aln",
        start_tree="output/strains/{strain}/besttree.nwk"
    output:
        final_tree="output/strains/{strain}/gubbins.final_tree.tre",
    group:
        "gubbins"
    log:
        "../../../logs/gubbins_{strain}.log"
    params:
        prefix=config['gubbins']['prefix'],
        tree_builder=config['gubbins']['tree_builder'],
        min_snp=config['gubbins']['min_snps'],
        min_window=config['gubbins']['min_window_size'],
        max_window=config['gubbins']['max_window_size'],
        iterations=config['gubbins']['iterations'],
    threads:
        2
    shell:
        "pushd output/strains/{wildcards.strain}/ && \
        run_gubbins.py map_variants.aln --prefix {params.prefix} \
        --starting-tree besttree.nwk --tree-builder {params.tree_builder} \
        --min-snps {params.min_snp} --min-window-size {params.min_window} \
        --max-window-size {params.max_window} --iterations {params.iterations} \
        --threads {threads} > {log} \
        && popd"

rule bactdating:
    input:
        final_tree="output/strains/{strain}/gubbins.final_tree.tre",
        metadata=config["transmission_metadata"]
    output:
        rds="output/strains/{strain}/bactdate_data.rds",
        tree="output/strains/{strain}/bactdate_tree.nwk",
        sorted="output/strains/{strain}/sorted.rds"
    group:
        "outbreak"
    params:
        script=config["poppipe_location"] + "/scripts/run_bactDating.R"
    log:
        "logs/bactdating_{strain}.log"
    shell:
        "Rscript --vanilla {params.script} output/strains/{wildcards.strain}/gubbins {input.metadata} {output.sorted} {output.tree} {output.rds} > {log}"

# dummy target rule
rule transmission:
    input:
        expand("output/strains/{strain}/transphylo_results.rds", strain=included_strain_ids)
    output:
        "done.txt"
    shell:
        "touch done.txt"

rule transphylo:
    input:
        rds="output/strains/{strain}/bactdate_data.rds",
        sorted="output/strains/{strain}/sorted.rds"
    output:
        rds="output/strains/{strain}/transphylo_results.rds"
    group:
        "outbreak"
    params:
        script=config["poppipe_location"] + "/scripts/run_transPhylo.R",
        w_shape=config['transphylo']['w_shape'],
        w_scale=config['transphylo']['w_scale'],
        mcmcIterations=config['transphylo']['mcmcIterations'],
        startNeg=config['transphylo']['startNeg'],
        startOff_r=config['transphylo']['startOff_r'],
        startOff_p=config['transphylo']['startOff_p'],
        startPi=config['transphylo']['startPi'],
        optiStart=config['transphylo']['optiStart'],
        dateT=config['transphylo']['dateT'],
        gubbins="output/strains/{strain}/gubbins"
    log:
        "logs/transphylo_{strain}.log"
    shell:
        "Rscript --vanilla {params.script} --rds {input.rds} --sorted {input.sorted} --output {output} \
        --gubbins {params.gubbins} --wshape {params.w_shape} \
        --wscale {params.w_scale} --mcmcIterations {params.mcmcIterations} \
        --startNeg {params.startNeg} --startOffr {params.startOff_r} \
        --startOffp {params.startOff_p} --startPi {params.startPi} \
        --optiStart {params.optiStart} --dateT {params.dateT} > {log}"

rule fastbaps:
    input:
        tree="output/strains/{strain}/besttree.nwk",
        align="output/strains/{strain}/align_variants.aln"
    output:
        "output/strains/{strain}/fastbaps_clusters.txt"
    group:
        "bapsclust"
    params:
        script=config["poppipe_location"] + "/scripts/run_fastbaps.R",
        levels=config['fastbaps']['levels']
    log:
        "logs/fastbaps_{strain}.log"
    threads:
        4
    shell:
        "Rscript --vanilla {params.script} {input.align} {output} {params.levels} {input.tree} {threads} > {log}"

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
    script:
        config["poppipe_location"] + "/scripts/tree_graft.py"

rule generate_dot:
    input:
        config["poppunk_h5"]
    output:
        dot="output/viz/mandrake.embedding.dot"
    group:
        "viz"
    params:
        db_prefix = db_prefix,
        n_samples = len(samples),
        perplexity = str(config['mandrake']['perplexity']),
        knn = str(config['mandrake']['knn']),
        maxIter = str(config['mandrake']['maxIter'])
    log:
        "logs/mandrake.log"
    threads:
        64
    script:
        config["poppipe_location"] + "/scripts/run_mandrake.py"

# use microreact api
rule make_microreact:
    input:
        tree="output/full_tree.nwk",
        clusters="output/all_clusters.txt",
        dot="output/viz/mandrake.embedding.dot"
    output:
        url="output/microreact_url.txt",
        microreact="output/json.microreact"
    params:
        example_file=config["poppipe_location"] + "/microreact_example.pkl",
        microreact_name=config['microreact']['name'],
        microreact_email=config['microreact']['email'],
        microreact_website=config['microreact']['website'],
        microreact_token=config['microreact']['api_token']
    group:
        "viz"
    script:
        config["poppipe_location"] + "/scripts/make_microreact.py"
