# Wrapper for Picard AddOrReplaceReadGroups

This wrapper provides a more structured interface for Picard's AddOrReplaceReadGroups tool by allowing read group parameters to be specified individually.

## Required Parameters

* `ID`: Read Group ID (required)
* `LB`: Read Group Library (required)
* `PL`: Read Group platform (e.g. ILLUMINA, SOLID) (required)
* `PU`: Read Group platform unit (required)
* `SM`: Read Group sample name (required)

## Optional Parameters

* `CN`: Read Group sequencing center name
* `DS`: Read Group description
* `DT`: Read Group run date
* `FO`: Read Group flow order
* `KS`: Read Group key sequence
* `PG`: Read Group program group
* `PI`: Read Group predicted insert size
* `PM`: Read Group platform model

## Additional Parameters

* `extra`: Any additional parameters to pass to Picard AddOrReplaceReadGroups

## Example

```python
rule add_readgroups:
    input:
        "mapped.bam"
    output:
        "mapped_with_rg.bam"
    params:
        # Required Parameters
        ID="flow_cell.lane",
        LB="library_prep_1",
        PL="ILLUMINA",
        PU="flow_cell.1.lane",
        SM="sample1",
        # Optional Parameters
        CN="sequencing_center",
        PM="NextSeq2000",
        # Additional Parameters
        extra="--VALIDATION_STRINGENCY LENIENT"
    wrapper:
        "file://path/to/damlab-wrappers/picard/addorreplacereadgroups"
```

## Author
* Will Dampier (@wnd) 