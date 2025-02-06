# Wrapper for Tree Rerooting

A wrapper that uses [DendroPy](https://dendropy.org/) to reroot phylogenetic trees at a specified taxon.

## Input
- Phylogenetic tree file in Newick format

## Output
- Rerooted tree file in Newick format

## Params
* `root_taxon` (required)
    The name of the taxon to root the tree at

## Example

```python
rule reroot_tree:
    input:
        "tree.nwk"
    output:
        "rerooted_tree.nwk"
    params:
        root_taxon="OutgroupSpecies"
    wrapper:
        "file:path/to/damlab-wrappers/phylo/reroot"

```

## Authors
* Will Dampier, PhD

## Software Requirements
* [DendroPy](https://dendropy.org/) (tested with v4.5.2)

## Notes
- The specified root taxon must exist in the tree
- The tree must be in Newick format 