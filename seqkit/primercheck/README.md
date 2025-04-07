# Seqkit Primercheck

This wrapper uses seqkit amplicon to check primers against sequence reads and produces a summary CSV file.
For each read, it reports the amplicon length for each primer pair (or None if no match).

## Input
* `reads`: Input sequence file (FASTA/Q) or BAM file
* `primers`: (optional) Tab-delimited primer file with columns:
    1. Primer name
    2. Forward primer (5'-3')
    3. Reverse primer (5'-3')

## Output
* Required:
    - CSV file with:
        - Rows: read IDs
        - Columns: primer names
        - Values: amplicon lengths or "None" if no match

* Optional:
    - Summary YAML file with:
        - Total number of sequences checked
        - Number of hits for each primer
        - Pairwise matrix of dual hits between primers

## Params
* `extra`: Optional parameters passed to seqkit amplicon
* `region`: Optional region string for BAM input (e.g. "chr1:1000-2000")
* `max_mismatch`: Maximum allowed mismatches in primer matching (default: 1)
* `primer_sets`: (optional) List or string of pre-defined primer sets to use.
  * `silicano-hiv` : 
    * Psi-F - CAGGACTCGGCTTGCTGAAG
    * Psi-R - GCACCCATCTCTCTCCTTCTAGC
    * Env-F - AGTGGTGCAGAGAGAAAAAAGAGC
    * Env-R - GTCTGGCCTGTACCGTCAGC
    * alt-Psi-F - GCAGGACTCGGCTTGCTG
    * alt-Psi-R - GCACCCATCTCTCTCTCCTTCTAG
  * `jones-hiv` :
    * Psi-F - CAGGACTCGGCTTGCTGAAG
    * Psi-R - GCACCCATCTCTCTCCTTCTAGC
    * Env-F - AGTGGTGCAGAGAGAAAAAAGAGC
    * Env-R - GTCTGGCCTGTACCGTCAGC
    * alt-Env-F - ACTATGGGCGCAGCGTC
    * alt-Env-R - CCCCAGACTGTGAGTTGCA
  * `deeks-hiv-subb` :
    * Psi-F - TCTCGACGCAGGACTCG
    * Psi-R - TACTGACGCTCTCGCACC
    * Env-F - AGTGGTGCAGAGAGAAAAAAGAGC
    * Env-R - GTCTGGCCTGTACCGTCAGC
  * `deeks-hiv-subc` :
    * Psi-F - TCTCGACGCAGGACTCG
    * Psi-R - TATTGACGCTCTCGCACC
    * Env-F - AGTGGTGGAGAGAGAAAAAAGAGC
    * Env-R - GTCTGGCCTGTACCGTCAGC

## Example
```python
# Basic usage with threads
rule check_primers:
    input:
        reads="reads.fastq",
        primers="primers.txt"
    output:
        "primer_results.csv"
    threads: 8  # Use 8 threads for processing
    wrapper:
        "file:path/to/damlab-wrappers/seqkit/primercheck"

# With BAM input and threads
rule check_primers_bam:
    input:
        reads="aligned.bam",
        primers="primers.txt"
    output:
        "primer_results.csv"
    params:
        region="chr1:1000-2000"
    threads: 8  # Use 8 threads for BAM processing and seqkit
    wrapper:
        "file:path/to/damlab-wrappers/seqkit/primercheck"

# With summary statistics
rule check_primers_with_summary:
    input:
        reads="reads.fastq",
        primers="primers.txt"
    output:
        "primer_results.csv",
        summary="primer_summary.yaml"
    wrapper:
        "file:path/to/damlab-wrappers/seqkit/primercheck"

# Basic usage with custom mismatch threshold
rule check_primers:
    input:
        reads="reads.fastq",
        primers="primers.txt"
    output:
        "primer_results.csv"
    params:
        max_mismatch=1  # More stringent matching
    wrapper:
        "file:path/to/damlab-wrappers/seqkit/primercheck"
```

## Authors
* Will Dampier, PhD

## Software Requirements
* [seqkit](https://bioinf.shenwei.me/seqkit/) (tested with v2.3.1)
* [samtools](http://www.htslib.org/) (for BAM input)
* Python packages:
    - pyyaml 