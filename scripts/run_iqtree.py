import subprocess
from shutil import copyfile
from pandas import read_table
import os

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

# change align_variants.aln
new_names = []
replacement_needed = False
for name in og_files:
    og_name=name.split()[0]
    if '#' in og_name:
        new_names.append(og_name.replace('#', '_'))
        replacement_needed = True
    else:
        new_names.append(og_name)

if replacement_needed:
    align_file=open(snakemake.input.alignment, 'r')
    new_file_name = snakemake.input.alignment.replace('align_variants.aln', 'temp_align_variants.aln')
    new_align_file = open(snakemake.params.temp, 'w')
    new_align_file = open(new_file_name, 'w')
    for line in align_file:
        line=line.strip()
        if ">" in line:
            new_align_file.write(">" + new_names.pop(0) + '\n')
        else:
            new_align_file.write(line + '\n')
    new_align_file.close()
    align_file.close()
    os.rename(new_file_name, snakemake.input.alignment)
