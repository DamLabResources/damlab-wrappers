# Phytreeviz Wrapper

A wrapper for [phytreeviz](https://github.com/damlab/phytreeviz), a tool for creating publication-quality visualizations of phylogenetic trees with customizable styling options.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Usage

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

```python
rule visualize_tree:
    input:
        "tree.nwk"
    output:
        "tree_viz.pdf"
    params:
        format="newick",           # Input tree format
        width=15,                  # Figure width in inches
        height=20,                 # Figure height in inches
        show_branch_support=True,  # Show bootstrap/support values
        color_scheme="viridis",    # Color scheme for branches
        dpi=300,                   # Resolution for raster formats
        extra=""                   # Additional parameters
    wrapper:
        "file:path/to/damlab-wrappers/phylo/phytreeviz"
```

## Parameters

- `format` (str, optional): Input tree format. Defaults to "newick".
- `width` (float, optional): Figure width in inches. Defaults to 10.
- `height` (float, optional): Figure height in inches. Defaults to 10.
- `show_branch_support` (bool, optional): Show bootstrap/support values. Defaults to False.
- `color_scheme` (str, optional): Color scheme for branches (e.g., "viridis", "magma"). Defaults to None.
- `dpi` (int, optional): Resolution for raster formats. Defaults to 300.
- `extra` (str, optional): Additional parameters to pass to phytreeviz. Defaults to "".
- `version` (str, optional): Specify a version of the wrapper to use.

## Input
* Phylogenetic tree file
  - Supports multiple formats (Newick, Nexus, PhyloXML)
  - Must be a valid tree file with proper formatting
  - Can include branch lengths and support values

## Output
* Tree visualization file
  - Supports multiple formats (PDF, PNG, SVG)
  - Publication-quality graphics
  - Customizable styling and layout

## Error Handling

The wrapper includes comprehensive error handling for:
- Missing input files
- Invalid parameter values
- Invalid tree formats
- Visualization errors

## Output Format
The output visualization includes:
- Properly scaled branches
- Bootstrap/support values (if enabled)
- Custom colors and styling
- High-resolution graphics

## Author
* Will Dampier, PhD

## Software Requirements
* [phytreeviz](https://github.com/damlab/phytreeviz) (tested with v0.2.0)

## License
This project is licensed under the MIT License - see the LICENSE file for details.