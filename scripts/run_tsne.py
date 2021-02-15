import subprocess
import os

if snakemake.config["tsne"]["use_gpu"]:
    subprocess.run("poppunk_sketch --query --ref-db " + snakemake.params['db_prefix'] + " --query-db " + snakemake.params['db_prefix'] + \
                    " --output " + snakemake.params['dist_prefix'] + " --read-k --use-gpu --gpu-id " + \
                    snakemake.params['device_id'] + "&> " + snakemake.log['sketch_log'],
                    shell=True, check=True)
else:
    subprocess.run("poppunk_sketch --query --ref-db " + snakemake.params['db_prefix'] + " --query-db " + snakemake.params['db_prefix'] + \
                    " --output " + snakemake.params['dist_prefix'] + " --read-k --cpus " + str(snakemake.threads) \
                     + " &> " + snakemake.log['sketch_log'],
                    shell=True, check=True)

subprocess.run("poppunk_tsne --distances " + snakemake.params['dist_prefix'] + " --output output --perplexity " + \
                snakemake.params['perplexity'] + " --verbosity 1 &> " + snakemake.log['tsne_log'],
                shell=True, check=True)
os.rename("output/output_perplexity" + snakemake.params['perplexity'] + "_accessory_tsne.dot",
          snakemake.output['tsne_out'])