import re
import subprocess

if snakemake.params["n_samples"] >= 100:
    subprocess.run(f"mandrake --sketches {snakemake.input} --output output/viz/mandrake" + \
               f" --use-accessory --perplexity {snakemake.params['perplexity']} --kNN {snakemake.params['knn']}" + \
               f" --maxIter {snakemake.params['maxIter']} --cpus {snakemake.threads} --no-clustering" + \
               f" 2> " + str(snakemake.log),
               shell=True, check=True)
else:
    for ofile in snakemake.output:
        subprocess.run(f"touch {ofile}", shell=True, check=True)
