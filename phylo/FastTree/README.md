# Wrapper for FastTree

[FastTree](http://www.microbesonline.org/fasttree/) infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences.

## Input
- Multiple sequence alignment file (FASTA format)

## Output
- Phylogenetic tree in Newick format

## Params
* `extra`
    Optional parameters that will be passed to FastTree (e.g., "-gtr -gamma")

## Example
```python
rule fasttree:
    input:
        "aligned.fasta"
    output:
        "tree.nwk"
    wrapper:
        "file:path/to/damlab-wrappers/phylo/FastTree"

rule fasttree:
    input:
        "aligned.fasta"
    output:
        "tree.nwk"
    params:
        extra="-gtr -gamma -boot 1000"
    wrapper:
        "file:path/to/damlab-wrappers/phylo/FastTree"
```

## Authors
* Will Dampier, PhD

## Software Requirements
* [FastTree](http://www.microbesonline.org/fasttree/) (tested with v2.1.11) 