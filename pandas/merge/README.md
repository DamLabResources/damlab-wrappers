# Pandas Merge Wrapper

This wrapper provides a convenient way to merge multiple CSV files using pandas merge functionality. It supports various merge types and join conditions, with proper handling of column suffixes for multiple file merges.

## Input
* List of CSV files to be merged

## Output
* Single merged CSV file

## Parameters
* `on` (optional)
    Column or index level names to join on. Can be:
    - A single column name (str) - used for all merges
    - A list of column names for all tables (List[str])
    - A list of lists where each inner list contains a key for that table (List[List[str]])
      In this case, the keys are used in a pairwise progressive pattern:
      For tables A, B, C with keys [['a'], ['b'], ['c']], the merges are:
      A.merge(B, left_on=['a'], right_on=['b'])
      result.merge(C, left_on=['b'], right_on=['c'])
* `how` (optional, default: "inner")
    Type of merge to be performed:
    - 'left': use only keys from left frame
    - 'right': use only keys from right frame
    - 'outer': use union of keys from both frames
    - 'inner': use intersection of keys from both frames
    - 'cross': creates the cartesian product
* `suffixes` (optional)
    Tuple of strings to append to overlapping column names.
    If not provided, defaults to ("_0", "_1", "_2", etc.) based on file order.
    **Important**: The number of suffixes must match the number of input files.

## Example
```python
rule merge_csv_files:
    input:
        csv_files=["file1.csv", "file2.csv", "file3.csv"]
    output:
        "merged_output.csv"
    params:
        on="common_column",
        how="outer",
        suffixes=("_first", "_second", "_third")  # Must provide one suffix per file
    wrapper:
        "file:path/to/damlab-wrappers/pandas/merge"
```

## Multiple File Merging
When merging multiple files, the wrapper requires one suffix per input file:
1. For N files, you must provide N suffixes
2. Each file will use its corresponding suffix when there are column name conflicts
3. If no suffixes are provided, the wrapper will automatically generate suffixes ("_0", "_1", etc.)

For example, with three files and suffixes=("_A", "_B", "_C"):
- Columns from file1 get suffix "_A"
- Columns from file2 get suffix "_B"
- Columns from file3 get suffix "_C"

## Different Column Names
When merging tables with different column names for their keys, you can use the pairwise progressive pattern:

```python
rule merge_different_columns:
    input:
        csv_files=["employees.csv", "salaries.csv", "positions.csv"]
    output:
        "merged_employees.csv"
    params:
        on=[["employee_id"], ["staff_id"], ["person_id"]],
        how="outer"
    wrapper:
        "file:path/to/damlab-wrappers/pandas/merge"
```

This will perform the following merges:
1. employees.csv.merge(salaries.csv, left_on=['employee_id'], right_on=['staff_id'])
2. result.merge(positions.csv, left_on=['staff_id'], right_on=['person_id'])

## Output Format
The output is a single CSV file containing the merged data from all input files, with appropriate suffixes on overlapping column names.

## Changelog
### v1.1.0
* Added support for pairwise progressive pattern in the `on` parameter
* Improved handling of different column names between tables
* Added better error handling and logging for merge operations

### v1.0.0
* Initial release with basic merge functionality
* Support for multiple file merging with suffixes
* Support for different merge types (inner, outer, left, right)

## Author
* Will Dampier

## Software Requirements
* [Pandas](https://pandas.pydata.org/) (tested with v2.0.0) 