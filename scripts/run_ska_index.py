import re
import subprocess

if re.match(r"\.(fq|fastq)\.$", str(snakemake.input['assembly'])):
    subprocess.run("ska fastq -o output/ska_index/" + str(snakemake.wildcards['sample']) + \
                   " -q " + str(snakemake.params['fastq_qual']) + " -c " + \
                   str(snakemake.params['fastq_cov']) + " " + str(snakemake.input['assembly']) + \
                   " > " + str(snakemake.log),
                   shell=True, check=True)
else:
    subprocess.run("ska fasta -o output/ska_index/" + str(snakemake.wildcards['sample']) + \
                   " " + str(snakemake.input['assembly']) + " > " + str(snakemake.log),
                   shell=True, check=True)