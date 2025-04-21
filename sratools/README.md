# SRA Toolkit Wrappers

This directory contains Snakemake wrappers for two SRA Toolkit commands:
1. `prefetch`: Downloads SRA files and their dependencies
2. `fasterq-dump`: Converts SRA files to FASTQ format

## Using the Wrappers Together

The recommended way to use these wrappers is in sequence, with the prefetch output feeding into fasterq-dump. Here's an example workflow:

```python
rule download_sra:
    input:
        # No input needed as we're downloading from SRA
    output:
        temp(directory("sra/SRR123456"))  # Mark as temporary to delete after processing
    params:
        accession="SRR123456",
        max_size="20G",
        transport="both",
        progress=True
    wrapper:
        "file:path/to/damlab-wrappers/sratools/prefetch"

rule convert_to_fastq:
    input:
        "sra/SRR123456"
    output:
        "fastq/SRR123456_1.fastq",
        "fastq/SRR123456_2.fastq"
    params:
        outdir="fastq",
        split_3=True,
        threads=8,
        progress=True
    wrapper:
        "file:path/to/damlab-wrappers/sratools/fasterq-dump"
```

For some SRA runs, this instances, this does not work and fails with an error of:
```bash
fasterq-dump.3.2.1 err: the input data is missing the QUALITY-column
```

In those instances, it is still possible to download the files using fasterq-dump directly as so:
```python
rule convert_to_fastq:
    input:
        # No input needed as we're using accession directly
    output:
        "test_output/DRR537798.fastq"
    params:
        accession="DRR537798",
        outdir="test_output",
        format='fastq',
        split_3=False,  # Don't split for single-ended data
        threads=4,
        progress=True
    wrapper:
        "file:path/to/damlab-wrappers/sratools/fasterq-dump"
```


## Important Considerations

### Disk Space Requirements

The SRA Toolkit requires significant disk space during processing:

1. **Prefetch Stage**:
   - Downloads the SRA file and its dependencies
   - Size varies by dataset but can be large (10s of GB)

2. **Fasterq-dump Stage**:
   - The final FASTQ files will be approximately 7 times the size of the SRA file
   - Temporary space (scratch space) of about 1.5 times the final FASTQ size is needed
   - Overall, you need approximately 17 times the size of the SRA file during conversion

Example calculation for a 10GB SRA file:
- Final FASTQ files: ~70GB
- Temporary space needed: ~105GB
- Total space needed during processing: ~175GB

### Using Temporary Files

It's highly recommended to mark the SRA directory as temporary using Snakemake's `temp()` function:
```python
output:
    temp(directory("sra/SRR123456"))
```

This ensures that:
1. The SRA files are automatically deleted after successful FASTQ conversion
2. Disk space is freed up for subsequent steps
3. The workflow remains clean and efficient

### Error Handling and Resumption

Both wrappers support error handling and resumption:

1. **Prefetch**:
   - Use `resume="yes"` to continue interrupted downloads
   - Use `verify="yes"` to ensure downloaded files are complete
   - If prefetch fails, running the same command again will resume the download

2. **Fasterq-dump**:
   - Use `force=True` to overwrite existing files if needed
   - The tool will automatically handle interruptions and can be restarted

## Best Practices

1. **Space Management**:
   - Always check available disk space before starting
   - Use `df -h` to monitor space during processing
   - Consider using a dedicated high-speed storage for temporary files

2. **Error Prevention**:
   - Set appropriate `max_size` in prefetch to prevent partial downloads
   - Use `verify="yes"` to ensure data integrity
   - Monitor system resources during processing

3. **Workflow Integration**:
   - Use Snakemake's `temp()` function for SRA files
   - Consider using `checkpoint` for large datasets
   - Implement proper error handling in your workflow

## Example Complete Workflow

```python
rule all:
    input:
        "fastq/SRR123456_1.fastq",
        "fastq/SRR123456_2.fastq"

rule download_sra:
    input:
        # No input needed as we're downloading from SRA
    output:
        temp(directory("sra/SRR123456"))
    params:
        accession="SRR123456",
        max_size="20G",
        transport="both",
        progress=True,
        verify="yes",
        resume="yes"
    wrapper:
        "file:path/to/damlab-wrappers/sratools/prefetch"

rule convert_to_fastq:
    input:
        "sra/SRR123456"
    output:
        "fastq/SRR123456_1.fastq",
        "fastq/SRR123456_2.fastq"
    params:
        outdir="fastq",
        split_3=True,
        threads=8,
        progress=True,
        temp="/scratch/temp"  # Use fast storage for temporary files
    wrapper:
        "file:path/to/damlab-wrappers/sratools/fasterq-dump"
```

## Additional Resources

- [SRA Toolkit Documentation](https://github.com/ncbi/sra-tools/wiki)
- [Snakemake Documentation](https://snakemake.readthedocs.io/)
- [SRA Toolkit HowTo Guide](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump) 