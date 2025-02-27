# Dorado Demultiplexing

Wrapper for the ONT dorado demultiplexing tool.
This tool splits reads into separate files based on the barcode sequence.

## Inputs
* `reads`: Path to input reads file or directory (required)

## Outputs
* `directory`: Output directory for demultiplexed reads

## Parameters
* `kit_name`: Name of the barcoding kit (required unless no_classify=True)
* `dorado_path`: Path to dorado executable (default: "dorado")
* `sample_sheet`: Path to sample sheet file (optional)
* `max_reads`: Maximum number of reads to process (optional)
* `read_ids`: Path to file containing read IDs to process (optional)
* `emit_fastq`: Output FASTQ instead of BAM (default: False)
* `emit_summary`: Generate summary file (default: False)
* `barcode_both_ends`: Require barcodes at both ends (default: False)
* `no_trim`: Disable barcode trimming (default: False)
* `sort_bam`: Sort output BAM files (default: False)
* `barcode_arrangement`: Custom barcode arrangement (optional)
* `barcode_sequences`: Custom barcode sequences (optional)
* `barcode_to_output`: Dictionary mapping barcodes to output filenames
* `tempdir`: Directory for temporary files (optional)
* `no_classify`: Skip barcode classification (default: False)

## Example

```python
rule dorado_demux:
    input:
        reads="path/to/reads.bam"
    output:
        directory("path/to/output_dir")
    params:
        kit_name="SQK-RBK004",
        emit_fastq=True,
        barcode_to_output={
            "barcode01": "sample1.fastq",
            "barcode02": "sample2.fastq"
        }
    threads: 8
    wrapper:
        "file://path/to/damlab-wrappers/dorado/demux" 
```

## Authors
* Will Dampier, PhD

## Software Requirements
* [DORADO](https://github.com/nanoporetech/dorado) (tested with v0.9.1)