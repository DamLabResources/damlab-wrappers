# FastTree Wrapper

A wrapper for [FastTree](http://www.microbesonline.org/fasttree/), a program for inferring approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Usage

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

```python
rule fasttree:
    input:
        "aligned.fasta"
    output:
        "tree.nwk"
    params:
        gtr=True,        # Use GTR model (for nucleotides)
        gamma=True,      # Use gamma-distributed rates
        boot=1000,       # Number of bootstrap replicates
        nt=True,         # Input is nucleotide sequences
        extra=""         # Additional parameters
    wrapper:
        "file:path/to/damlab-wrappers/phylo/FastTree"
```

## Parameters

- `gtr` (bool, optional): Use generalized time-reversible model (for nucleotides). Defaults to False.
- `gamma` (bool, optional): Use gamma-distributed rates across sites. Defaults to False.
- `boot` (int, optional): Number of bootstrap replicates. Defaults to None (no bootstrapping).
- `nt` (bool, optional): Input is nucleotide sequences. Defaults to True.
- `extra` (str, optional): Additional parameters to pass to FastTree. Defaults to "".
- `version` (str, optional): Specify a version of the wrapper to use.

## Input
* Multiple sequence alignment file in FASTA format
  - Can be nucleotide or protein sequences
  - Must be properly aligned (all sequences same length)

## Output
* Phylogenetic tree in Newick format
  - Contains branch lengths
  - Contains bootstrap values (if bootstrapping enabled)
  - Can be visualized with tools like FigTree or iTOL

## Error Handling

The wrapper includes comprehensive error handling for:
- Missing input files
- Invalid parameter values
- FastTree execution errors

## Output Format
The output is a Newick format tree file, which includes:
- Branch lengths representing evolutionary distances
- Bootstrap support values (if bootstrapping enabled)
- Internal node labels (if any)

## Author
* Will Dampier, PhD

## Software Requirements
* [FastTree](http://www.microbesonline.org/fasttree/) (tested with v2.1.11)

## License
This project is licensed under the MIT License - see the LICENSE file for details. 