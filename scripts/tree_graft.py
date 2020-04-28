import dendropy

njtree = dendropy.Tree.get(path=snakemake.input.nj_tree, 
                           schema="newick", 
                           preserve_underscores=True)

for mltree_file in snakemake.input['ml_trees']:
    mltree = dendropy.Tree.get(path=mltree_file, 
                               schema="newick", 
                               preserve_underscores=True)
    
    for sample in mltree.taxon_namespace:
        top_node = njtree.find_node_with_taxon_label(sample.label)
        if top_node:
            matched_node = mltree.find_node_with_taxon_label(sample.label)
            mltree.reroot_at_edge(matched_node.edge, update_bipartitions=False)
            top_node.add_child(mltree.nodes()[0])
            break

#njtree.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
#njtree.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
njtree.write(path=str(snakemake.output),
            schema="newick",
            suppress_rooting=True,
            unquoted_underscores=True)