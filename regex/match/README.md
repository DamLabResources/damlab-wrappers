# Regex Pattern Matching Wrapper

This wrapper searches for regex patterns in FASTA files and outputs the matches to a CSV file.

## Input

- FASTX (optionally .gz) file containing sequences to search
- List of regex patterns to search for

## Output

- CSV file containing:
  - Read name
  - One column per pattern containing the first match found (if any)
- Optional metrics file in YAML format

## Parameters

- `patterns`: List of regex patterns to search for in the sequences
  - Each pattern will be compiled with BESTMATCH and IGNORECASE flags
  - The first match for each pattern will be recorded
  - **Note**: When using regex patterns with curly braces (e.g., `{5}`) in Snakemake, you must:
    1. Double the curly braces (`{{5}}`) to escape them from Snakemake's wildcard expansion
    2. Wrap the patterns in a lambda function to prevent early evaluation
    ```python
    params:
        patterns=lambda wildcards: [
            r"ATG[ATGC]{{5}}TAG",  # Correct: double curly braces
            r"GAT[ATGC]{{5}}CTA"   # Correct: double curly braces
        ]
    ```

## Example Usage

```python
rule regex_match:
    input:
        fasta="input.fasta"
    output:
        csv="matches.csv",
        metrics="metrics.yaml"
    params:
        patterns=lambda wildcards: [
            r"ATG[ATGC]{{10,20}}TAG",  # Example pattern 1
            r"GAT[ATGC]{{5,15}}CTA"    # Example pattern 2
        ]
    wrapper:
        "damlab-wrappers/regex/match"
```

## Output Format

The CSV output will have the following columns:
- `read_name`: The name of the read from the FASTA file
- `pattern_1`: First match for the first pattern (if found)
- `pattern_2`: First match for the second pattern (if found)
- ... and so on for each pattern

## Metrics

If a metrics file is specified, it will contain:
- Total number of reads processed
- Number of matches found for each pattern 