# Wrapper for phytreeviz

[phytreeviz](https://github.com/damlab/phytreeviz) is a tool for visualizing phylogenetic trees.

## Input
- Phylogenetic tree file (typically in Newick format)

## Output
- Tree visualization (PDF, PNG, or other supported formats)

## Params
* `extra`
    Optional parameters that will be passed to phytreeviz
* `format`
    Input tree format (default: "newick")
* `width`
    Figure width in inches (default: 10)
* `height`
    Figure height in inches (default: 10)

## Example
```python
rule visualize_tree:
    input:
        "tree.nwk"
    output:
        "tree_viz.pdf"
    wrapper:
        "file:path/to/damlab-wrappers/phylo/phytreeviz"


rule visualize_tree:
    input:
        "tree.nwk"
    output:
        "tree_viz.pdf"
    params:
        format="newick",
        width=15,
        height=20,
        extra="--show_branch_support --color_scheme viridis"
    wrapper:
        "file:path/to/damlab-wrappers/phylo/phytreeviz"

```

## Authors
* Will Dampier, PhD

## Software Requirements
* [phytreeviz](https://github.com/damlab/phytreeviz) (tested with v0.2.0)