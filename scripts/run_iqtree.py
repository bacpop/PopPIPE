import subprocess
from shutil import copyfile

from ete3 import Tree

def midpoint_root(infile, outfile):
    t = Tree(infile)
    t.set_outgroup(t.get_midpoint_outgroup())
    t.write(format=5, outfile=outfile) # format 5: internal and leaf branches + leaf names

if snakemake.params.enabled:
    iqtree_cmd = "iqtree --quiet -s " + snakemake.input.alignment + " -t " + snakemake.input.start_tree + \
                 " -T " + str(snakemake.threads) + " --prefix " + snakemake.params.prefix
    if snakemake.params.mode == "full":
        iqtree_cmd += " -m " + snakemake.params.model
    elif snakemake.params.mode == "fast":
        iqtree_cmd += " --fast"

    subprocess.run(iqtree_cmd, shell=True, check=True)
else:
    copyfile(snakemake.input.start_tree, snakemake.output.unrooted)

midpoint_root(snakemake.output.unrooted, snakemake.output.rooted)
