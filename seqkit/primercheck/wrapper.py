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
from typing import Dict, Set, Optional, Union, List, Tuple
import pysam

if "snakemake" not in locals():
    #  Keeps linters happy but doesn't impact function
    import snakemake # type: ignore

# Predefined primer sets
PRIMER_SETS = {
    'silicano-hiv': {
        'Psi-F': 'CAGGACTCGGCTTGCTGAAG',
        'Psi-R': 'GCACCCATCTCTCTCCTTCTAGC',
        'Env-F': 'AGTGGTGCAGAGAGAAAAAAGAGC',
        'Env-R': 'GTCTGGCCTGTACCGTCAGC',
        'alt-Psi-F': 'GCAGGACTCGGCTTGCTG',
        'alt-Psi-R': 'GCACCCATCTCTCTCTCCTTCTAG',
    },
    'jones-hiv': {
        'Psi-F': 'CAGGACTCGGCTTGCTGAAG',
        'Psi-R': 'GCACCCATCTCTCTCCTTCTAGC',
        'Env-F': 'AGTGGTGCAGAGAGAAAAAAGAGC',
        'Env-R': 'GTCTGGCCTGTACCGTCAGC',
        'alt-Env-F': 'ACTATGGGCGCAGCGTC',
        'alt-Env-R': 'CCCCAGACTGTGAGTTGCA'
    },
    'deeks-hiv-subb': {
        'Psi-F': 'TCTCGACGCAGGACTCG',
        'Psi-R': 'TACTGACGCTCTCGCACC',
        'Env-F': 'AGTGGTGCAGAGAGAAAAAAGAGC',
        'Env-R': 'GTCTGGCCTGTACCGTCAGC'
    },
    'deeks-hiv-subc': {
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

def parse_amplicon_results(bed_file: Path, all_primers: Optional[Set[str]] = None) -> tuple[Dict, Set]:
    """Parse seqkit amplicon BED output into results dictionary.
    
    Args:
        bed_file: Path to BED format output from seqkit amplicon
        all_primers: Optional set of all primer names that should be included
    
    Returns:
        Tuple of (results dict, set of primer names)
    """
    results = defaultdict(dict)
    found_primers = set()
    
    with open(bed_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            read_id = fields[0]
            primer_name = fields[3]
            if len(fields) > 6:
                amplicon_length = len(fields[6])
            else:
                amplicon_length = None
            
            results[read_id][primer_name] = amplicon_length
            found_primers.add(primer_name)
    
    # Include all expected primers in the set
    primer_names = found_primers
    if all_primers:
        primer_names = all_primers
    
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

def get_all_primer_names(primer_file: Optional[Path], primer_sets: Optional[Union[str, List[str]]]) -> Set[str]:
    """Get set of all primer names from input file and/or primer sets.
    
    Args:
        primer_file: Optional path to input primer file
        primer_sets: Optional string or list of primer set names
    
    Returns:
        Set of all primer names
    """
    primer_names = set()
    
    # Get primers from file
    if primer_file:
        with open(primer_file) as f:
            for line in f:
                name = line.strip().split('\t')[0]
                primer_names.add(name)
    
    # Get primers from predefined sets
    if primer_sets:
        if isinstance(primer_sets, str):
            primer_sets = [primer_sets]
        
        for set_name in primer_sets:
            if set_name not in PRIMER_SETS:
                raise ValueError(f"Unknown primer set: {set_name}")
            
            primers_dict = PRIMER_SETS[set_name]
            # Group primers into pairs (F/R)
            for name in primers_dict:
                base_name = name[:-2]  # Remove -F or -R
                if name.endswith('-F'):
                    primer_names.add(f"{set_name}_{base_name}")
    
    return primer_names

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
    
    # Count hits per primer (include zeros for unused primers)
    primer_hits = {primer: sum(1 for read in results.values() if primer in read) 
                  for primer in primer_names}
    
    # Generate pairwise matrix (include zeros for all combinations)
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

def is_bam_file(filename: str) -> bool:
    """Check if a file is a BAM file by looking at the first few bytes.
    
    Args:
        filename: Path to the file to check
        
    Returns:
        True if file appears to be BAM format
    """
    return filename.endswith('.bam')

def parse_region(region: str) -> Tuple[str, int, int]:
    """Parse a region string into reference name, start, and end positions.
    
    Args:
        region: Region string in format "chr:start-end"
        
    Returns:
        Tuple of (reference name, start position, end position)
        
    Raises:
        ValueError: If region string is not properly formatted
    """
    try:
        chrom, pos = region.split(':')
        start, end = map(int, pos.split('-'))
        return chrom, start, end
    except ValueError:
        raise ValueError(f"Invalid region format: {region}. Expected format: chr:start-end")

def is_read_in_region(read: pysam.AlignedSegment, ref_name: str, start: int, end: int) -> bool:
    """Check if a read overlaps with the specified region.
    
    Args:
        read: Pysam AlignedSegment object
        ref_name: Reference sequence name
        start: Region start position (1-based)
        end: Region end position (1-based)
        
    Returns:
        True if read overlaps with region
    """
    # Skip unmapped reads
    if read.is_unmapped:
        return False
    
    # Check if read is on the correct reference
    if read.reference_name != ref_name:
        return False
    
    # Convert to 0-based coordinates for comparison
    read_start = read.reference_start  # Already 0-based
    read_end = read.reference_end or read_start  # Handle reads with no end position
    
    # Convert region to 0-based
    region_start = start - 1
    region_end = end
    
    # Check for overlap
    return read_start < region_end and read_end > region_start

def extract_reads_from_bam(
    bam_file: Path,
    output_fasta: Path,
    threads: int = 1,
    region: Optional[str] = None,
    log: str = ""
) -> None:
    """Extract reads from BAM file to FASTA format.
    
    Args:
        bam_file: Path to input BAM file
        output_fasta: Path to write FASTA output
        threads: Number of threads to use
        region: Optional region string (e.g. "chr1:1000-2000")
        log: Log redirect string for shell command
    """
    # Parse region if provided
    region_info = None
    if region:
        region_info = parse_region(region)
    
    # Open BAM file
    with pysam.AlignmentFile(str(bam_file), "rb", threads=threads) as bam:
        # Open output FASTA
        with open(output_fasta, 'w') as fasta:
            # Iterate through all reads
            for read in bam:
                # Skip secondary/supplementary alignments
                if read.is_secondary or read.is_supplementary:
                    continue
                
                # Check region if specified
                if region_info and not is_read_in_region(read, *region_info):
                    continue
                
                # Get sequence
                seq = read.get_forward_sequence()
                
                # Write to FASTA
                fasta.write(f">{read.query_name}\n{seq}\n")

# Extract arguments from snakemake object
reads = snakemake.input.reads
primers = snakemake.input.get("primers", None)
output_csv = snakemake.output[0]
summary_yaml = snakemake.output.get("summary", None)

# Get optional parameters
extra = snakemake.params.get("extra", "")
sample_name = snakemake.params.get("sample_name", None)
primer_sets = snakemake.params.get("primer_sets", None)
region = snakemake.params.get("region", None)
max_mismatch = snakemake.params.get("max_mismatch", 1)  # Default to 1 mismatch

# Get threads
threads = snakemake.threads

# Setup log
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

with tempfile.TemporaryDirectory() as temp_dir:
    temp_primers = Path(temp_dir) / "primers.txt"
    temp_bed = Path(temp_dir) / "amplicons.bed"
    
    # Handle BAM input if necessary
    input_reads = reads
    if is_bam_file(reads):
        temp_fasta = Path(temp_dir) / "reads.fasta"
        extract_reads_from_bam(Path(reads), temp_fasta, threads, region, log)
        input_reads = temp_fasta
    
    # Create primer file and get all expected primer names
    create_primer_file(temp_primers, primers, primer_sets)
    all_primers = get_all_primer_names(primers, primer_sets)
    
    # Run seqkit amplicon with max_mismatch parameter
    shell(
        f"seqkit amplicon --bed -m {max_mismatch} -j {threads} --primer-file {temp_primers} {input_reads} > {temp_bed} {log}"
    )
    
    # Parse results with all expected primers
    results, primer_names = parse_amplicon_results(temp_bed, all_primers)
    
    # Write CSV output
    write_csv_results(Path(output_csv), results, primer_names)
    
    # Generate and write summary if requested
    if summary_yaml:
        summary = generate_summary(results, primer_names, sample_name)
        with open(summary_yaml, 'w') as f:
            f.write("# Seqkit Primer Check Summary\n")
            yaml.dump(summary, f, default_flow_style=False)
