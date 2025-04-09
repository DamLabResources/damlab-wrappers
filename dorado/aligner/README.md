# Dorado Aligner

Wrapper for the ONT dorado aligner tool, which uses minimap2 to align basecalled reads.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Installation

Dorado can be installed from the Nanopore GitHub repository:

```bash
# Clone the repository
git clone https://github.com/nanoporetech/dorado.git
cd dorado

# Build dorado (requires CMake and a C++ compiler)
mkdir build && cd build
cmake ..
make -j

# Add dorado to your PATH
export PATH=$PATH:$(pwd)/bin
```

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

## Usage

```python
rule dorado_align:
    input:
        index="path/to/mm2.mmi",
        calls="path/to/basecalls.bam"
    output:
        "path/to/aligned.bam"
    params:
        k=15,
        x="map-ont",
        secondary=True,
        mm2_opts="-A2 -B4",
        emit_summary=True,
        bed_file="annotations.bed"
    threads: 8
    wrapper:
        "file://path/to/damlab-wrappers/dorado/aligner"
```

## Parameters

- `k` (int, optional): Minimap2 k-mer size for alignment (max 28)
- `w` (int, optional): Minimap2 minimizer window size
- `I` (int, optional): Minimap2 index batch size
- `N` (int, optional): Maximum number of secondary alignments to retain
- `r` (int, optional): Chaining/alignment bandwidth
- `x` (str, optional): Minimap2 preset (default: "lr:hq")
- `secondary` (bool, optional): Output secondary alignments
- `Y` (bool, optional): Use soft clipping for supplementary alignments
- `junc_bed` (str, optional): Gene annotations in BED12 format
- `mm2_opts` (str, optional): Additional minimap2 options
- `dorado_path` (str, optional): Path to dorado executable (default: "dorado")
- `recursive` (bool, optional): Search subfolders for input files
- `output_dir` (str, optional): Output directory for files
- `emit_summary` (bool, optional): Emit alignment summary file
- `bed_file` (str, optional): BED file for overlap counting
- `max_reads` (int, optional): Maximum number of reads to process
- `verbose` (bool, optional): Increase verbosity

## Input
* `index`: Path to minimap2 index file (required)
* `calls`: Path to basecalled reads in BAM format (required)

## Output
* BAM file containing aligned reads

## Error Handling

The wrapper includes error handling for:
- Missing input files
- Invalid parameters
- Dorado execution errors
- Alignment failures

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Authors
* Will Dampier, PhD

## Software Requirements
* [DORADO](https://github.com/nanoporetech/dorado) (tested with v0.9.1) 