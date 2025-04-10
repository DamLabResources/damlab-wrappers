# Seqkit Primercheck Wrapper

A wrapper for the [seqkit amplicon](https://bioinf.shenwei.me/seqkit/usage/#amplicon) tool that checks primers against sequence reads and produces a summary CSV file. For each read, it reports the amplicon length for each primer pair (or None if no match).

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Usage

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

```python
# Basic usage with threads
rule check_primers:
    input:
        reads="reads.fastq",
        primers="primers.txt"
    output:
        "primer_results.csv"
    params:
        version="1.0.0"  # Wrapper version
    threads: 8  # Use 8 threads for processing
    wrapper:
        "file://path/to/damlab-wrappers/seqkit/primercheck"

# With BAM input and region specification
rule check_primers_bam:
    input:
        reads="aligned.bam",
        primers="primers.txt"
    output:
        "primer_results.csv"
    params:
        region="chr1:1000-2000",
        version="1.0.0"
    threads: 8  # Use 8 threads for BAM processing and seqkit
    wrapper:
        "file://path/to/damlab-wrappers/seqkit/primercheck"

# With summary statistics
rule check_primers_with_summary:
    input:
        reads="reads.fastq",
        primers="primers.txt"
    output:
        "primer_results.csv",
        summary="primer_summary.yaml"
    params:
        version="1.0.0"
    wrapper:
        "file://path/to/damlab-wrappers/seqkit/primercheck"

# Using predefined primer sets
rule check_primers_with_sets:
    input:
        reads="reads.fastq"
    output:
        "primer_results.csv"
    params:
        primer_sets=["silicano-hiv", "jones-hiv"],
        version="1.0.0"
    wrapper:
        "file://path/to/damlab-wrappers/seqkit/primercheck"
```

## Parameters

### Input Options
- `reads` (str): Input sequence file (FASTA/Q) or BAM file
- `primers` (str, optional): Tab-delimited primer file with columns:
  1. Primer name
  2. Forward primer (5'-3')
  3. Reverse primer (5'-3')

### Output
- `output.csv` (str): CSV file with:
  - Rows: read IDs
  - Columns: primer names
  - Values: amplicon lengths or "None" if no match
- `summary` (str, optional): YAML file with:
  - Total number of sequences checked
  - Number of hits for each primer
  - Pairwise matrix of dual hits between primers

### Optional Parameters
- `extra` (str, optional): Additional parameters to pass to seqkit amplicon
- `region` (str, optional): Region string for BAM input (e.g. "chr1:1000-2000")
- `max_mismatch` (int, optional): Maximum allowed mismatches in primer matching (default: 1)
- `primer_sets` (str or list, optional): Pre-defined primer sets to use
- `sample_name` (str, optional): Name of the sample for summary statistics
- `version` (str, optional): Specify a version of the wrapper to use

## Predefined Primer Sets

The wrapper includes several predefined primer sets for HIV analysis:

### silicano-hiv
- Psi-F: CAGGACTCGGCTTGCTGAAG
- Psi-R: GCACCCATCTCTCTCCTTCTAGC
- Env-F: AGTGGTGCAGAGAGAAAAAAGAGC
- Env-R: GTCTGGCCTGTACCGTCAGC
- alt-Psi-F: GCAGGACTCGGCTTGCTG
- alt-Psi-R: GCACCCATCTCTCTCTCCTTCTAG

### jones-hiv
- Psi-F: CAGGACTCGGCTTGCTGAAG
- Psi-R: GCACCCATCTCTCTCCTTCTAGC
- Env-F: AGTGGTGCAGAGAGAAAAAAGAGC
- Env-R: GTCTGGCCTGTACCGTCAGC
- alt-Env-F: ACTATGGGCGCAGCGTC
- alt-Env-R: CCCCAGACTGTGAGTTGCA

### deeks-hiv-subb
- Psi-F: TCTCGACGCAGGACTCG
- Psi-R: TACTGACGCTCTCGCACC
- Env-F: AGTGGTGCAGAGAGAAAAAAGAGC
- Env-R: GTCTGGCCTGTACCGTCAGC

### deeks-hiv-subc
- Psi-F: TCTCGACGCAGGACTCG
- Psi-R: TATTGACGCTCTCGCACC
- Env-F: AGTGGTGGAGAGAGAAAAAAGAGC
- Env-R: GTCTGGCCTGTACCGTCAGC

## Input
* Sequence reads (FASTA/Q or BAM format)
  - Must be valid sequence files
  - BAM files are automatically converted to FASTA for processing
* Primer file (optional if using predefined sets)
  - Tab-delimited format with primer name, forward sequence, and reverse sequence
  - Each line represents a primer pair

## Output
* CSV file with primer check results
  - One row per read
  - One column per primer pair
  - Values indicate amplicon length or "None" if no match
* Optional YAML summary file
  - Total sequence count
  - Hit counts per primer
  - Pairwise hit matrix

## Error Handling

The wrapper includes comprehensive error handling for:
- Missing input files
- Invalid primer files
- Unknown primer sets
- Invalid region specifications
- BAM file processing errors
- Command execution errors

## Notes
- The wrapper supports both FASTA/Q and BAM input formats
- BAM files are automatically converted to FASTA for processing
- Region specification is only applicable to BAM input
- Predefined primer sets can be used instead of a primer file
- Multiple primer sets can be combined

## Author
* Will Dampier, PhD

## Software Requirements
* [seqkit](https://bioinf.shenwei.me/seqkit/) (tested with v2.3.1)
* [samtools](http://www.htslib.org/) (for BAM input)
* Python packages:
  - pyyaml
  - pysam
  - Python >= 3.10

## License
This project is licensed under the MIT License - see the LICENSE file for details. 