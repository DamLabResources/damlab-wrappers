# Tree Rerooting Wrapper

A wrapper for [DendroPy](https://dendropy.org/)'s tree rerooting functionality, allowing users to reroot phylogenetic trees at a specified taxon. This is particularly useful for ensuring proper tree orientation and outgroup placement.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Usage

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

```python
rule reroot_tree:
    input:
        "tree.nwk"
    output:
        "rerooted_tree.nwk"
    params:
        root_taxon="OutgroupSpecies",      # Required: taxon to root at
        schema="newick",                    # Input/output tree format
        preserve_branch_lengths=True,       # Keep branch lengths
        preserve_support_values=True,       # Keep bootstrap values
        version="1.0.0"                     # Wrapper version to use
    wrapper:
        "file:path/to/damlab-wrappers/phylo/reroot"
```

## Parameters

- `root_taxon` (str, required): The name of the taxon to root the tree at.
- `schema` (str, optional): Tree file format. Defaults to "newick".
- `preserve_branch_lengths` (bool, optional): Whether to keep branch lengths. Defaults to True.
- `preserve_support_values` (bool, optional): Whether to keep bootstrap/support values. Defaults to True.
- `version` (str, optional): Specify a version of the wrapper to use.

## Input
* Phylogenetic tree file
  - Supports multiple formats (primarily Newick)
  - Must contain the specified root taxon
  - Can include branch lengths and support values

## Output
* Rerooted tree file
  - Same format as input
  - Preserves branch lengths (if requested)
  - Preserves support values (if requested)
  - Properly oriented with specified root taxon

## Error Handling

The wrapper includes comprehensive error handling for:
- Missing input files
- Missing required parameters
- Root taxon not found in tree
- Invalid tree formats
- Tree reading/writing errors
- Rerooting errors

## Notes
- The specified root taxon must exist in the tree
- Underscores in taxon names are preserved
- Branch lengths and support values can be optionally preserved
- Tree bipartitions are updated during rerooting

## Author
* Will Dampier, PhD

## Software Requirements
* [DendroPy](https://dendropy.org/) (tested with v4.5.2)

## License
This project is licensed under the MIT License - see the LICENSE file for details. 