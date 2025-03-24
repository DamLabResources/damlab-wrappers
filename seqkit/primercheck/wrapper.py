"""Wrapper for seqkit amplicon to check primers"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

from snakemake.shell import shell # type: ignore
import csv
from pathlib import Path
from collections import defaultdict
import yaml
from itertools import combinations
import tempfile
from typing import Dict, Set, Optional, Union, List

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact function
    import snakemake # type: ignore

# Predefined primer sets
PRIMER_SETS = {
    'silicano-hiv': {
        'Psi-F': 'CAGGACTCGGCTTGCTGAAG',
        'Psi-R': 'GCTAGAAGGAGAGAGATGGGTGC',
        'Env-F': 'AGTGGTGCAGAGAGAAAAAAGAGC',
        'Env-R': 'GTCTGGCCTGTACCGTCAGC',
        'alt-Psi-F': 'GCAGGACTCGGCTTGCTG',
        'alt-Psi-R': 'CTAGAAGGAGAGAGAGATGGGTGC'
    },
    'jones-hiv': {
        'Psi-F': 'CAGGACTCGGCTTGCTGAAG',
        'Psi-R': 'GCACCCATCTCTCTCCTTCTAGC',
        'Env-F': 'AGTGGTGCAGAGAGAAAAAAGAGC',
        'Env-R': 'GTCTGGCCTGTACCGTCAGC',
        'alt-Env-F': 'ACTATGGGCGCAGCGTC',
        'alt-Env-R': 'CCCCAGACTGTGAGTTGCA'
    },
    'Deeks-hiv-SubB': {
        'Psi-F': 'TCTCGACGCAGGACTCG',
        'Psi-R': 'TACTGACGCTCTCGCACC',
        'Env-F': 'AGTGGTGCAGAGAGAAAAAAGAGC',
        'Env-R': 'GTCTGGCCTGTACCGTCAGC'
    },
    'Deeks-hiv-SubC': {
        'Psi-F': 'TCTCGACGCAGGACTCG',
        'Psi-R': 'TATTGACGCTCTCGCACC',
        'Env-F': 'AGTGGTGGAGAGAGAAAAAAGAGC',
        'Env-R': 'GTCTGGCCTGTACCGTCAGC'
    }
}

def create_primer_file(
    output_path: Path,
    primer_file: Optional[Path] = None,
    primer_sets: Optional[Union[str, List[str]]] = None
) -> None:
    """Create a combined primer file from input file and/or preset primer sets.
    
    Args:
        output_path: Path to write the combined primer file
        primer_file: Optional path to input primer file
        primer_sets: Optional string or list of primer set names to include
    
    Raises:
        ValueError: If neither primer_file nor primer_sets is provided
        ValueError: If specified primer set doesn't exist
    """
    if primer_file is None and primer_sets is None:
        raise ValueError("Either primer_file or primer_sets must be provided")
    
    with open(output_path, 'w') as f:
        # Add primers from file if provided
        if primer_file:
            with open(primer_file) as pf:
                for line in pf:
                    f.write(line)
        
        # Add primers from predefined sets
        if primer_sets:
            # Convert string to list if single set provided
            if isinstance(primer_sets, str):
                primer_sets = [primer_sets]
            
            for set_name in primer_sets:
                if set_name not in PRIMER_SETS:
                    raise ValueError(f"Unknown primer set: {set_name}")
                
                primers_dict = PRIMER_SETS[set_name]
                primer_pairs = {}
                
                # Group primers into pairs (F/R)
                for name, seq in primers_dict.items():
                    base_name = name[:-2]  # Remove -F or -R
                    suffix = name[-2:]  # Get -F or -R
                    if suffix in ('-F', '-R'):
                        if base_name not in primer_pairs:
                            primer_pairs[base_name] = {}
                        primer_pairs[base_name][suffix] = seq
                
                # Write complete pairs to file
                for base_name, pair in primer_pairs.items():
                    if '-F' in pair and '-R' in pair:
                        f.write(f"{set_name}_{base_name}\t{pair['-F']}\t{pair['-R']}\n")

def parse_amplicon_results(bed_file: Path) -> tuple[Dict, Set]:
    """Parse seqkit amplicon BED output into results dictionary.
    
    Args:
        bed_file: Path to BED format output from seqkit amplicon
    
    Returns:
        Tuple of (results dict, set of primer names)
    """
    results = defaultdict(dict)
    primer_names = set()
    
    with open(bed_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            read_id = fields[0]
            primer_name = fields[3]
            amplicon_length = len(fields[6])
            
            results[read_id][primer_name] = amplicon_length
            primer_names.add(primer_name)
    
    return results, primer_names

def write_csv_results(
    output_path: Path,
    results: Dict,
    primer_names: Set
) -> None:
    """Write amplicon results to CSV file.
    
    Args:
        output_path: Path to output CSV file
        results: Dictionary of results from parse_amplicon_results
        primer_names: Set of primer names from parse_amplicon_results
    """
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header
        header = ['read_id'] + sorted(list(primer_names))
        writer.writerow(header)
        
        # Write data
        for read_id in sorted(results.keys()):
            row = [read_id]
            for primer in header[1:]:
                row.append(results[read_id].get(primer, 'None'))
            writer.writerow(row)

def generate_summary(
    results: Dict,
    primer_names: Set,
    sample_name: Optional[str] = None
) -> Dict:
    """Generate summary statistics from amplicon results.
    
    Args:
        results: Dictionary of results from parse_amplicon_results
        primer_names: Set of primer names from parse_amplicon_results
        sample_name: Optional name of the sample
    
    Returns:
        Dictionary of summary statistics
    """
    total_seqs = len(results)
    
    # Count hits per primer
    primer_hits = {primer: sum(1 for read in results.values() if primer in read) 
                  for primer in primer_names}
    
    # Generate pairwise matrix
    pairwise_hits = {}
    for p1, p2 in combinations(sorted(primer_names), 2):
        dual_hits = sum(1 for read in results.values() 
                       if p1 in read and p2 in read)
        pairwise_hits[f"{p1}_x_{p2}"] = dual_hits
    
    summary = {
        'total_sequences': total_seqs,
        'primer_hits': primer_hits,
        'pairwise_hits': pairwise_hits
    }
    
    if sample_name:
        summary['sample_name'] = sample_name
    
    return summary


# Extract arguments from snakemake object
reads = snakemake.input.reads
primers = snakemake.input.get("primers", None)
output_csv = snakemake.output[0]
summary_yaml = snakemake.output.get("summary", None)

# Get optional parameters
extra = snakemake.params.get("extra", "")
sample_name = snakemake.params.get("sample_name", None)
primer_sets = snakemake.params.get("primer_sets", None)

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

with tempfile.TemporaryDirectory() as temp_dir:
    temp_primers = Path(temp_dir) / "primers.txt"
    temp_bed = Path(temp_dir) / "amplicons.bed"
    
    # Create primer file
    create_primer_file(temp_primers, primers, primer_sets)
    
    # Run seqkit amplicon
    shell(
        f"seqkit amplicon --bed --primer-file {temp_primers} {reads} > {temp_bed} {log}"
    )
    
    # Parse results
    results, primer_names = parse_amplicon_results(temp_bed)
    
    # Write CSV output
    write_csv_results(Path(output_csv), results, primer_names)
    
    # Generate and write summary if requested
    if summary_yaml:
        summary = generate_summary(results, primer_names, sample_name)
        with open(summary_yaml, 'w') as f:
            f.write("# Seqkit Primer Check Summary\n")
            yaml.dump(summary, f, default_flow_style=False)
