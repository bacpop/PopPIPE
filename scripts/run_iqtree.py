import subprocess
from shutil import copyfile

if snakemake.params.enabled:
    iqtree_cmd = "iqtree -s " + snakemake.input.alignment + " -t " snakemake.input.start_tree +
                 " -T " + snakemake.threads + " --prefix " + snakemake.output.unrooted
    if snakemake.params.mode == "full":
        iqtree_cmd += " -m " + snakemake.params.model
    elif snakemake.params.mode == "fast": 
        iqtree_cmd += " --fast"
    iqtree_cmd += " > " + snakemake.log

    subprocess.run(iqtree_cmd, shell=True, check=True)
else:
    copyfile(snakemake.input.start_tree, output.unrooted)