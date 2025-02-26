# Dorado Aligner

Wrapper for the ONT dorado aligner tool, which uses minimap2 to align basecalled reads.

## Inputs
* `index`: Path to minimap2 index file (required)
* `calls`: Path to basecalled reads in BAM format (required)

## Outputs
* BAM file containing aligned reads

## Parameters
* `k`: Minimap2 k-mer size for alignment (max 28)
* `w`: Minimap2 minimizer window size
* `I`: Minimap2 index batch size
* `N`: Maximum number of secondary alignments to retain
* `r`: Chaining/alignment bandwidth
* `x`: Minimap2 preset (default: "lr:hq")
* `secondary`: Output secondary alignments (boolean)
* `Y`: Use soft clipping for supplementary alignments (boolean)
* `junc_bed`: Gene annotations in BED12 format
* `mm2_opts`: Additional minimap2 options
* `dorado_path`: Path to dorado executable (default: "dorado")
* `recursive`: Search subfolders for input files (boolean)
* `output_dir`: Output directory for files
* `emit_summary`: Emit alignment summary file (boolean)
* `bed_file`: BED file for overlap counting
* `max_reads`: Maximum number of reads to process
* `verbose`: Increase verbosity (boolean)

## Example

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

## Authors
* Will Dampier, PhD

## Software Requirements
* [DORADO](https://github.com/nanoporetech/dorado) (tested with v0.9.1) 