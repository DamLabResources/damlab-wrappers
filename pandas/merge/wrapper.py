"""Wrapper for pandas merge functionality"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import pandas as pd
import yaml
import logging
from typing import List, Optional
from pathlib import Path

# Configure basic logging to stdout
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('pandas-merge-wrapper')

if "snakemake" not in locals():
    import snakemake  # type: ignore

def merge_csv_files(input_files: List[str], 
                   on: Optional[str | List[str]] = None,
                   how: str = 'inner',
                   suffixes: Optional[tuple] = None) -> pd.DataFrame:
    """Merge multiple CSV files iteratively.
    
    Args:
        input_files: List of CSV file paths to merge
        on: Column(s) to merge on
        how: Type of merge to perform ('left', 'right', 'outer', 'inner')
        suffixes: Suffixes to apply to overlapping column names
        
    Returns:
        Merged DataFrame
    """
    if not input_files:
        raise ValueError("No input files provided")
        
    logger.info(f"Merging {len(input_files)} CSV files")
    logger.info(f"Merge parameters - on: {on}, how: {how}")
    
    # Read first file as base
    result = pd.read_csv(input_files[0])
    logger.info(f"Read base file: {input_files[0]} with shape {result.shape}")
    
    # If suffixes not provided, create default suffixes based on number of files
    if suffixes is None:
        suffixes = tuple(f"_{i}" for i in range(len(input_files)))
    
    # Enforce that the number of suffixes matches the number of files
    assert len(suffixes) == len(input_files), f"Number of suffixes ({len(suffixes)}) must match number of input files ({len(input_files)})"
    
    # Iteratively merge remaining files
    for i, file in enumerate(input_files[1:], 1):
        df = pd.read_csv(file)
        logger.info(f"Merging file: {file} with shape {df.shape}")
        
        result = result.merge(
            df,
            on=on,
            how=how,
            suffixes=(suffixes[i-1], suffixes[i])
        )
        logger.info(f"Shape after merge: {result.shape}")
    
    return result

def generate_metrics(df: pd.DataFrame, input_files: List[str], params: dict) -> dict:
    """Generate metrics about the merge operation.
    
    Args:
        df: Merged DataFrame
        input_files: List of input file paths
        params: Dictionary of merge parameters
        
    Returns:
        Dictionary of metrics
    """
    metrics = {
        'input_files': input_files,
        'parameters': params,
        'output_shape': {
            'rows': df.shape[0],
            'columns': df.shape[1]
        },
        'column_names': list(df.columns),
        'memory_usage_bytes': df.memory_usage(deep=True).sum()
    }
    
    return metrics

def main():
    # Get input files and output file
    input_files = snakemake.input.csv_files
    output_file = snakemake.output[0]
    metrics_file = snakemake.output.get('metrics', None)  # Get optional metrics file
    
    # Get parameters with defaults
    on = snakemake.params.get("on", None)
    how = snakemake.params.get("how", "inner")
    suffixes = snakemake.params.get("suffixes", None)
    
    params = {
        'on': on,
        'how': how,
        'suffixes': suffixes
    }
    
    try:
        # Perform the merge
        merged_df = merge_csv_files(
            input_files,
            on=on,
            how=how,
            suffixes=suffixes
        )
        
        # Save merged result
        merged_df.to_csv(output_file, index=False)
        logger.info(f"Saved merged result to {output_file}")
        
        # Generate and save metrics if requested
        if metrics_file:
            metrics = generate_metrics(merged_df, input_files, params)
            
            with open(metrics_file, 'w') as f:
                f.write("# CSV Merge Metrics\n")
                yaml.dump(metrics, f, default_flow_style=False)
            logger.info(f"Saved metrics to {metrics_file}")
            
    except Exception as e:
        logger.error(f"Error merging files: {e}")
        raise

if __name__ == "__main__":
    main() 