import subprocess
from shutil import copyfile
from pandas import read_table
import os

from ete3 import Tree

def midpoint_root(infile, outfile):
    t = Tree(infile)
    t.set_outgroup(t.get_midpoint_outgroup())
    t.write(format=5, outfile=outfile) # format 5: internal and leaf branches + leaf names

def iqtree_cmd(alignment, start_tree, threads, prefix, model, mode):
    cmd = f"iqtree --quiet -st DNA -s {alignment} -t {start_tree} -T {threads} --prefix {prefix}"
    if mode == "full":
        cmd += f" -m {model}"
    elif mode == "fast":
        cmd += " --fast"
    return cmd

if snakemake.params.enabled:
    cmd = iqtree_cmd(snakemake.input.alignment, snakemake.input.start_tree, str(snakemake.threads), snakemake.params.prefix, snakemake.params.model, snakemake.params.mode)
    ret_code = subprocess.run(cmd, shell=True)
    # Automatically retry with varsites if ASC fails
    if ret_code.returncode != 0:
        cmd = iqtree_cmd(snakemake.params.prefix + ".varsites.phy", snakemake.input.start_tree, str(snakemake.threads), snakemake.params.prefix, snakemake.params.model, snakemake.params.mode)
        subprocess.run(cmd, shell=True, check=True)
else:
    copyfile(snakemake.input.start_tree, snakemake.output.unrooted)

midpoint_root(snakemake.output.unrooted, snakemake.output.rooted)

# Change any hashes in names back from underscores
rooted_file=open(snakemake.output.rooted, 'r')
samples = read_table(snakemake.input.rfiles, header=None, sep='\t')
og_files = samples.iloc[:, 0].values.tolist()

for line in rooted_file:
    for name in og_files:
        og_name=name.split()[0]
        if '#' in og_name:
            new_name = og_name.replace('#', '_')
            line = line.replace(new_name, og_name)
rooted_file.close()

file = open(snakemake.output.rooted, 'w')
file.write(line)
file.close()
