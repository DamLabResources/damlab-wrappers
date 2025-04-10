# Picard AddOrReplaceReadGroups Wrapper

A wrapper for [Picard](https://broadinstitute.github.io/picard/)'s AddOrReplaceReadGroups tool, providing a structured interface for adding or replacing read groups in BAM files.

## Version

Current version: 1.0.0 (First stable release)

For a detailed list of changes, see the [CHANGELOG.md](CHANGELOG.md).

## Usage

For general information on how to use wrappers in Snakemake, please refer to the [root README.md](../../../README.md).

```python
rule add_readgroups:
    input:
        "mapped.bam"
    output:
        "mapped_with_rg.bam"
    params:
        # Required Parameters
        ID="flow_cell.lane",        # Read Group ID
        LB="library_prep_1",        # Read Group Library
        PL="ILLUMINA",              # Read Group platform
        PU="flow_cell.1.lane",      # Read Group platform unit
        SM="sample1",               # Read Group sample name
        
        # Optional Parameters
        CN="sequencing_center",     # Read Group sequencing center
        PM="NextSeq2000",           # Read Group platform model
        
        # Additional Parameters
        validation_stringency="LENIENT",  # Validation stringency
        create_index=True,                # Create BAM index
        compression_level=5,              # BAM compression level
        extra="",                         # Additional parameters
        version="1.0.0"                   # Wrapper version
    wrapper:
        "file://path/to/damlab-wrappers/picard/addorreplacereadgroups"
```

## Parameters

### Required Parameters
- `ID` (str): Read Group ID (required)
- `LB` (str): Read Group Library (required)
- `PL` (str): Read Group platform (e.g. ILLUMINA, SOLID) (required)
- `PU` (str): Read Group platform unit (required)
- `SM` (str): Read Group sample name (required)

### Optional Parameters
- `CN` (str, optional): Read Group sequencing center name
- `DS` (str, optional): Read Group description
- `DT` (str, optional): Read Group run date
- `FO` (str, optional): Read Group flow order
- `KS` (str, optional): Read Group key sequence
- `PG` (str, optional): Read Group program group
- `PI` (int, optional): Read Group predicted insert size
- `PM` (str, optional): Read Group platform model
- `validation_stringency` (str, optional): Validation stringency (SILENT, LENIENT, STRICT). Defaults to None.
- `create_index` (bool, optional): Create BAM index. Defaults to True.
- `compression_level` (int, optional): BAM compression level (0-9). Defaults to None.
- `extra` (str, optional): Additional parameters to pass to Picard. Defaults to "".
- `version` (str, optional): Specify a version of the wrapper to use.

## Input
* BAM file
  - Must be a valid BAM file
  - Can be sorted or unsorted
  - Can have existing read groups or not

## Output
* BAM file with read groups added/replaced
  - Contains all reads from input BAM
  - Read groups updated according to parameters
  - BAM index created (if create_index=True)

## Error Handling

The wrapper includes comprehensive error handling for:
- Missing input files
- Missing required parameters
- Invalid parameter values
- Picard execution errors

## Notes
- Read group parameters are critical for downstream analysis
- The ID parameter must be unique for each read group
- The SM parameter should match the sample name used in variant calling
- Validation stringency can be adjusted if input BAM has issues

## Author
* Will Dampier, PhD

## Software Requirements
* [Picard](https://broadinstitute.github.io/picard/) (tested with v2.27.5)
* Java runtime environment

## License
This project is licensed under the MIT License - see the LICENSE file for details. 