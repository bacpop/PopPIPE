import dendropy

def midpoint_root(infile, outfile):
    from ete3 import Tree
    t = Tree(infile)
    t.set_outgroup(t.get_midpoint_outgroup())
    t.write(format=5, outfile=outfile) # format 5: internal and leaf branches + leaf names

fulltree = dendropy.Tree.get(path=snakemake.input.overall_tree,
                           schema="newick",
                           preserve_underscores=True)

for (mltree_file, njtree_file) in zip(snakemake.input['ml_trees'], snakemake.input['nj_trees']):
    mltree = dendropy.Tree.get(path=mltree_file,
                               schema="newick",
                               preserve_underscores=True)
    njtree = dendropy.Tree.get(path=njtree_file,
                               schema="newick",
                               preserve_underscores=True)
    mltree.scale_edges(njtree.length() / fulltree.length())

    for sample in mltree.taxon_namespace:
        top_node = fulltree.find_node_with_taxon_label(sample.label)
        if top_node:
            matched_node = mltree.find_node_with_taxon_label(sample.label)
            mltree.reroot_at_edge(matched_node.edge, update_bipartitions=False)

            # Set a non null edge length for join
            mltree.nodes()[0]._set_edge_length(0)
            mltree.nodes()[1]._set_edge_length(0)
            top_node.add_child(mltree.nodes()[0])
            break

# Get rid of taxa at internal nodes
for node in fulltree.internal_nodes():
    node.taxon = None

# As tree struct has been messed with not all functions work
# Hack is to go to/from file which sorts this out
fulltree.write(path=str(snakemake.output),
            schema="newick",
            suppress_rooting=True,
            unquoted_underscores=True)

midpoint_root(str(snakemake.output), str(snakemake.output))
