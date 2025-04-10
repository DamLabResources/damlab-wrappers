# ClipQS Wrapper

A Snakemake wrapper for clipping and orienting sequences using minimap2.
This tool aligns sequences to a reference, clips unaligned regions, and ensures proper orientation.

## Version

Current version: 1.0.0

See [CHANGELOG.md](CHANGELOG.md) for detailed changes.

## Usage

```python
rule clip_sequences:
    input:
        sequences = "path/to/sequences.fasta",
        reference = "path/to/reference.fasta"
    output:
        "path/to/processed.fasta"
    params:
        min_coverage = 0.2,  # Optional: Minimum coverage threshold (default: 0.2)
        include_reference = True,  # Optional: Include reference in output (default: True)
        version = "1.0.0"  # Optional: Wrapper version to use
    log:
        "logs/clipqs/{sample}.log"
    wrapper:
        "file://path/to/damlab-wrappers/strainline/clipqs"
```

## Parameters

### Required Parameters

None - all parameters have sensible defaults.

### Optional Parameters

- `min_coverage`: Minimum coverage threshold for keeping sequences. Must be between 0 and 1. Default: 0.2
- `include_reference`: Whether to include the reference sequence in the output. Default: True
- `version`: Wrapper version to use. Default: "1.0.0"

## Input

- `sequences`: FASTA file containing sequences to process
- `reference`: FASTA file containing the reference sequence for alignment

## Output

- FASTA file containing processed sequences with:
  - Clipped unaligned regions
  - Proper orientation matching reference
  - Coverage information in sequence headers
  - Optional inclusion of reference sequence

## Error Handling

The wrapper includes comprehensive error handling for common issues:

- Missing input files
- Invalid parameter values
- Empty reference file
- Failed aligner initialization
- Sequence processing errors
- Output file writing errors

## Notes

- The wrapper uses minimap2's "map-ont" preset, optimized for long-read sequences
- Sequences with coverage below the threshold are filtered out
- Coverage information is added to sequence headers
- The reference sequence can be optionally included in the output

## Author

Will Dampier (wnd22@drexel.edu)

## Software Requirements

- [minimap2](https://github.com/lh3/minimap2)
- [mappy](https://pypi.org/project/mappy/) Python package

## License

MIT 