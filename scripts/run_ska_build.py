import re
import subprocess

ss = "--single-strand" if snakemake.params['single_strand'] else ""

subprocess.run(f"ska build -o {snakemake.params['skf_prefix']}" + \
               f" -f {snakemake.input['samples']} -k {snakemake.params['kmer']}" + \
               f" --min-qual {snakemake.params['fastq_qual']} --min-count {snakemake.params['fastq_cov']}" + \
               f" {ss} --threads {snakemake.threads}" + \
               f" -v 2> " + str(snakemake.log),
               shell=True, check=True)

