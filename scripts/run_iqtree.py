import subprocess
from shutil import copyfile
from pandas import *
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
#
# # Change any hashes in names back from underscores
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
# -----------------------------------------------------------------------------
# # change align_variants.aln
# align_file=open(snakemake.input.alignment, 'r')
# new_align_file = open(snakemake.params.temp, 'w')
# for line in align_file:
#     line=line.strip()
#     if ">" in line:
#         for name in og_files:
#             og_name=name.split()[0]
#             if '#' in og_name:
#                 new_name = og_name.replace('#', '_')
#                 line = line.replace(new_name, og_name)
#         new_align_file.write(line + '\n')
#         # new_align_file.write('\n')
#     else:
#         new_align_file.write(line + '\n')
#         # new_align_file.write('\n')
# new_align_file.close()
#
# # os.remove(snakemake.input.alignment)
# os.rename(snakemake.params.temp ,snakemake.input.alignment)
# # delete align_variants
# # rename temp_align_variants to align_variants
