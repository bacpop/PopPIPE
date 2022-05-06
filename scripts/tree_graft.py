from ete3 import Tree

def midpoint_root(infile, outfile):
    t = Tree(infile)
    t.set_outgroup(t.get_midpoint_outgroup())
    t.write(format=5, outfile=outfile) # format 5: internal and leaf branches + leaf names

fulltree = Tree(snakemake.input.overall_tree)

for (mltree_file, njtree_file) in zip(snakemake.input['ml_trees'], snakemake.input['nj_trees']):
    mltree = Tree(path=mltree_file)
    njtree = Tree(path=njtree_file)
    ml_length = 0
    nj_length = 0
    for node in mltree.traverse():
        ml_length += node.dist
    for node in njtree.traverse():
        nj_length += node.dist
    scale = nj_length / ml_length
    for node in mltree.traverse():
        node.dist /= scale

    for sample in mltree.iter_leaves():
        top_node = fulltree.search_nodes(name=sample.name)
        if top_node:
            matched_node = mltree.find_node_with_taxon_label(sample.label)
            mltree.set_outgroup(matched_node)

            # Set a non null edge length for join
            top_node.add_child(mltree.nodes, dist=0)
            break

fulltree.set_outgroup(fulltree.get_midpoint_outgroup())
fulltree.write(format=5, outfile=str(snakemake.output))

