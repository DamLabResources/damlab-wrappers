# BAM2GEL

A wrapper for creating gel-like visualizations from BAM files using cigarmath. This tool processes one or more BAM files and creates a stacked histogram visualization that resembles an agarose gel, with each BAM file represented as a different lane.

## Features

- Process multiple BAM files in parallel
- Three different analysis modes:
  - `read_length`: Visualize read length distribution
  - `reference_block_size`: Visualize reference block size distribution
  - `query_block_size`: Visualize query block size distribution
- Memory-efficient streaming operations for large BAM files
- Customizable lane names and visualization parameters
- Outputs both visualization and detailed metrics

## Parameters

- `mode`: Analysis mode (default: "read_length")
  - Options: "read_length", "reference_block_size", "query_block_size"
- `names`: List of names for each lane (default: Sample_1, Sample_2, etc.)
- `bin_width`: Width of histogram bins in base pairs (default: 10)
- `max_size`: Maximum size to include in visualization (default: 1000)

## Outputs

1. PNG file containing the gel visualization
2. YAML file containing detailed metrics for each lane:
   - Total reads
   - Mean size
   - Median size
   - Minimum size
   - Maximum size

## Example Usage

```python
rule bam2gel:
    input:
        bams=["sample1.bam", "sample2.bam"]
    output:
        png="gel_visualization.png",
        yaml="gel_metrics.yaml"
    params:
        mode="read_length",
        names=["Control", "Treatment"],
        bin_width=20,
        max_size=2000
    wrapper:
        "damlab-wrappers/cigarmath/bam2gel"
```

## Dependencies

- cigarmath
- pandas
- seaborn
- matplotlib
- numpy
- pyyaml 