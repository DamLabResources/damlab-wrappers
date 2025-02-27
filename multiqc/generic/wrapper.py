"""Wrapper for creating generic MultiQC-compatible placeholder logs"""

import yaml
from snakemake.shell import shell
from pathlib import Path

if "snakemake" not in locals():
    import snakemake # type: ignore

def write_generic_log(params, output_file):
    """Write parameters in MultiQC-compatible YAML format."""
    with open(output_file, 'w') as handle:
        handle.write(f"# Generic MultiQC Log\n")
        # Convert all values to strings to ensure YAML compatibility
        formatted_params = {k: str(v) for k, v in params.items()}
        yaml.dump(formatted_params, handle, default_flow_style=False)

# Get output file
output_log = snakemake.output[0]

# Get all parameters from snakemake
params_dict = dict(snakemake.params)

# Add sample name if provided, otherwise use output filename stem
sample_name = params_dict.pop('sample_name', Path(output_log).stem)
params_dict['sample_name'] = sample_name

# Write the log file
write_generic_log(params_dict, output_log) 