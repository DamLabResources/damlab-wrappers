# Strainline Wrapper

A Snakemake wrapper for the Strainline haplotype reconstruction tool. Strainline is designed for strain-aware long-read haplotype phasing, supporting both ONT and PacBio platforms.

## Version

Current version: 1.0.0

See [CHANGELOG.md](CHANGELOG.md) for detailed changes.

## Usage

```python
rule strainline:
    input:
        reads = "path/to/reads.fastq"
    output:
        # Option 1: Output directory containing all results
        directory = directory("path/to/output_dir")
        # OR Option 2: Output single FASTA file with final haplotypes
        haplotypes = "path/to/haplotypes.fa"
    params:
        prefix = "/path/to/strainline/installation",  # Required: Path to Strainline installation
        platform = "ont",  # Optional: Sequencing platform ("ont" or "pb", default: "ont")
        extra_params = "",  # Optional: Additional Strainline parameters
        threads = 1,  # Optional: Number of threads to use
        version = "1.0.0"  # Optional: Wrapper version to use
    threads: 1
    log:
        "logs/strainline/{sample}.log"
    wrapper:
        "file://path/to/damlab-wrappers/seqkit/primercheck"

```

## Parameters

### Required Parameters

- `prefix`: Path to the Strainline installation directory. Must contain the Strainline script and Daccord installation.

### Optional Parameters

- `platform`: Sequencing platform to use. Must be either "ont" (Oxford Nanopore) or "pb" (PacBio). Default: "ont"
- `extra_params`: Additional parameters to pass to Strainline. Default: ""
- `threads`: Number of threads to use. Default: 1
- `version`: Wrapper version to use. Default: "1.0.0"

## Input

- `reads`: Input long-read sequences in FASTQ format.

## Output

The wrapper supports two output modes:

1. **Directory Mode**: Outputs a directory containing all Strainline results, including intermediate files and the final haplotypes.
2. **Haplotypes Mode**: Outputs a single FASTA file containing only the final haplotypes.

## Error Handling

The wrapper includes comprehensive error handling for common issues:

- Missing input files
- Invalid platform specification
- Missing or invalid Strainline installation
- Missing Daccord installation
- Missing haplotype output file after execution

## Notes

- The wrapper requires a complete Strainline installation at the specified prefix path
- The Daccord tool must be installed in the expected location within the Strainline installation
- For large datasets, consider increasing the number of threads for better performance
- The wrapper supports both ONT and PacBio data, but parameters may need to be adjusted based on the platform

## Author

Will Dampier (wnd22@drexel.edu)

## Software Requirements

- Strainline
- Daccord
- Snakemake

## License

MIT 