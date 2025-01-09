# Wrapper for jaspar2logo

A small utility for creating sequence logos from jaspar formatted files.
These are often produced by `grna/summarize_hits` but can be count tables from anywhere.

Currently, only a this is only a very simple wrapper with a lot of assumption.
More to come as needed.


# Wrapper application example:
```python
rule :
    input:
        jaspar = 'path/to/jaspar'
    output:
        logo = 'path/to/logo.png'
    params:
        dpi = 100,
        name_filter = None, # Specify a single name a multi-file
        # TODO # Specify target sequence and mark mismatches
        target_sequence = None, 
    wrapper:
        "http://10.11.19.24:3080/Data-Garden/snakemake-wrappers/raw/main/bio/jaspar2heatmap"
```


# Further documentation 
can be found [here](https://link.to.tool)
