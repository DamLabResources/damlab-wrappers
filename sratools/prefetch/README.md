# SRA Toolkit Prefetch Wrapper

Version: 1.0.0

This wrapper provides a Snakemake interface for the SRA Toolkit's `prefetch` command, which downloads SRA files and their dependencies.

## Input
* SRA accession number (e.g., SRR123456)

## Output
* Directory containing downloaded SRA files and dependencies

## Parameters
* `accession` (required)
    SRA accession number to download
* `type` (optional, default: "sra")
    File type to download
* `transport` (optional, default: "both")
    Transport method: "http", "fasp", or "both"
* `min_size` (optional)
    Minimum file size to download in KB
* `max_size` (optional, default: "20G")
    Maximum file size to download in KB
* `force` (optional, default: "no")
    Force object download: "no", "yes", "all", or "ALL"
* `resume` (optional, default: "yes")
    Resume partial downloads: "yes" or "no"
* `verify` (optional, default: "yes")
    Verify after download: "yes" or "no"
* `progress` (optional, default: true)
    Show progress
* `heartbeat` (optional, default: 1)
    Time period in minutes to display download progress
* `extra` (optional, default: "")
    Additional parameters to pass to prefetch

## Example
```python
rule download_sra:
    input:
        # No input needed as we're downloading from SRA
    output:
        directory("sra/SRR123456")
    params:
        accession="SRR123456",
        max_size="20G",
        transport="both",
        progress=True
    wrapper:
        "file:path/to/damlab-wrappers/sratools/prefetch"
```

## Output Format
The output is a directory containing the downloaded SRA files and their dependencies. The directory name will match the accession number.

## Author
* William Dam

## Software Requirements
* [SRA Toolkit](https://github.com/ncbi/sra-tools) (tested with v3.0.6) 