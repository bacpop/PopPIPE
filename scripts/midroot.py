import dendropy

tree = dendropy.Tree.get(path=snakemake.output.unrooted, schema="newick")
tree.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
tree.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
tree.write(path=str(snakemake.output.rooted),
            schema="newick",
            suppress_rooting=True,
            unquoted_underscores=True)
