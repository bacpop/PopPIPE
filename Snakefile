import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.1.2")

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

# table with sample, sequence, cluster
samples = pd.read_table(config["poppunk_rfile"]).set_index("sample", drop=False)
clusters = pd.read_table(config["poppunk_clusters"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")
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
        "{cluster}"
    output:
        "clusters/{cluster}/rfile.txt")
    script:
        "{config[script_loc]}/update_cluster_dirs.py" #TODO script to read PopPUNK results, create cluster dirs, add any new samples

# Use sketchlib to extract distances
# TODO write wrapper script which generates(/updates) input from a PopPUNK run
rule sketchlib_dists:
    input:
        database = config.poppunk_db,
        rfile = "clusters/{cluster}/rfiles.txt"
    output:
        "clusters/{cluster}/dists.npy",
        "clusters/{cluster}/dists.pkl"
    shell:
        "poppunk_sketch --query --ref-db {input.database} --query-db {input.database} --subset {input.rfile} "
        "--output clusters/{cluster}/dists --min-k {config[min_k]} --max-k {config[max_k]} --k-step {config[k_step]} "
        "--sketch-size {config[sketch_size]} --cpus {threads}"

# rapidnj
rule generate_nj:
    input:
        "clusters/{cluster}/dists.npy",
        "clusters/{cluster}/dists.pkl" 
    output:
        "clusters/{cluster}/njtree.nwk"
    script:
        "{config[script_loc]}/run_rapidnj.py"  #TODO use function from PopPUNK

# ska for alignment
rule ska_index:
    input:
        assembly=lambda wildcards: samples.loc[wildcards.sample, "sample_fa"]
    output:
        "clusters/{cluster}/{sample}.ska"
    shell:
        "ska fasta -o {output} {input.assembly}"

rule ska_align:
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

rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"


rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"


rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
