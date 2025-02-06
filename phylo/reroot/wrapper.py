"""Wrapper for tree rerooting using dendropy"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import dendropy # type: ignore

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact funtion
    import snakemake # type: ignore

# Extract arguments from snakemake object
input_tree = snakemake.input[0]
output_tree = snakemake.output[0]

# Get required root taxon parameter
root_taxon = snakemake.params.get("root_taxon")
if not root_taxon:
    raise ValueError("root_taxon parameter must be specified")

# Read the tree
tree = dendropy.Tree.get(
    path=input_tree,
    schema="newick"
)

# Find the node corresponding to the root taxon
root_node = None
for node in tree.leaf_node_iter():
    if node.taxon and node.taxon.label == root_taxon:
        root_node = node
        break

if not root_node:
    raise ValueError(f"Root taxon '{root_taxon}' not found in tree")

# Reroot the tree
tree.reroot_at_node(root_node)

# Write the rerooted tree
tree.write(
    path=output_tree,
    schema="newick"
) 