import subprocess
import dendropy
from shutil import copyfile

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

tree = dendropy.Tree.get(path=snakemake.output.unrooted, schema="newick")
tree.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
tree.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
tree.write(path=str(snakemake.output.rooted),
            schema="newick",
            suppress_rooting=True,
            unquoted_underscores=True)