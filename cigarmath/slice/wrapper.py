"""Wrapper for slicing aligned reads that overlap a specified region"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import pysam # type: ignore
import cigarmath as cm # type: ignore
import yaml
import re
import os
from typing import Tuple, Optional, Iterator

if "snakemake" not in locals():
    import snakemake # type: ignore

def parse_region(region_str: str) -> Tuple[str, int, int]:
    """Parse region string in format 'chr:start-end' or 'chr'."""
    match = re.match(r"(\w+)(?::(\d+)-(\d+))?", region_str)
    if not match:
        raise ValueError(f"Invalid region format: {region_str}. Expected 'chr:start-end' or 'chr'")
    
    chrom = match.group(1)
    start = int(match.group(2)) if match.group(2) else None
    end = int(match.group(3)) if match.group(3) else None
    
    return chrom, start, end

def check_bam_index(bam_path: str) -> bool:
    """Check if BAM index (.bai) file exists."""
    index_path = f"{bam_path}.bai"
    alt_index_path = f"{os.path.splitext(bam_path)[0]}.bai"
    return os.path.exists(index_path) or os.path.exists(alt_index_path)

def get_segment_stream(bam_path: str, chrom: Optional[str], min_mapq: int) -> Iterator:
    """Get a stream of aligned segments, using fetch if index exists."""
    has_index = check_bam_index(bam_path)
    
    if has_index and chrom:
        # Use fetch with chromosome if index exists
        return cm.io.segment_stream_pysam(
            bam_path, 
            mode='rb',
            fetch=chrom,  # Fetch just the chromosome to improve performance
            min_mapq=min_mapq
        )
    else:
        # No index or no specific chromosome, iterate through all reads
        return cm.io.segment_stream_pysam(
            bam_path, 
            mode='rb',
            min_mapq=min_mapq
        )

# Get input/output files
input_bam = snakemake.input[0]
output_fastq = snakemake.output[0]
output_metrics = snakemake.output[1] if len(snakemake.output) > 1 else None

# Get required region parameter
region = snakemake.params.region
chrom, region_start, region_end = parse_region(region)

# Get optional parameters
region_name = snakemake.params.get("region_name", region)
sample_name = snakemake.params.get("sample_name", "")
min_mapq = snakemake.params.get("min_mapq", 0)

# Initialize counters for metrics
metrics = {
    "region": region,
    "region_name": region_name,
    "sample_name": sample_name,
    "total_segments_processed": 0,
    "segments_overlapping_region": 0,
    "used_index": check_bam_index(input_bam)
}

# Stream segments from BAM file
segment_stream = get_segment_stream(input_bam, chrom, min_mapq)

# Open output FASTQ file
with open(output_fastq, 'w') as fastq_out:
    # Process each segment
    for segment in segment_stream:
        metrics["total_segments_processed"] += 1
        
        # Skip segments that don't map to our target chromosome
        if chrom and segment.reference_name != chrom:
            continue
            
        if not segment.cigartuples or not segment.query_sequence:
            continue
            
        # Check if segment overlaps the region
        segment_start = segment.reference_start
        segment_end = segment_start + cm.reference_offset(segment.cigartuples)
        
        # Skip if the segment doesn't overlap the region
        if region_end and segment_start >= region_end:
            continue
        if region_start and segment_end <= region_start:
            continue
            
        
        
        # Slice the segment to the region using cigar_iterator_reference_slice
        
        sliced_iterator = cm.cigar_iterator_reference_slice(
            segment.cigartuples,
            reference_start=segment.reference_start,
            region_reference_start=region_start,
            region_reference_end=region_end
        )
        
        #print(segment.query_sequence)
        # Attach query sequence information
        sliced_iterator = cm.iterator_attach(
            sliced_iterator,
            query_sequence=segment.query_sequence if segment.query_sequence else None,
            query_qualities=segment.query_qualities if segment.query_qualities else None
        )
        
        # Extract the sequence and quality from the sliced iterator
        seq_bases = []
        qual_scores = []
        
        try:
            for cigar_index in sliced_iterator:
                if cigar_index.query_letter is not None:
                    seq_bases.append(cigar_index.query_letter)
                if cigar_index.query_quality is not None:
                    qual_scores.append(cigar_index.query_quality)
        except IndexError:
            continue
            print(e)
            print('Query sequence length', len(segment.query_sequence))
            print('Query cigar length', cm.inferred_query_sequence_length(segment.cigartuples))

            print(segment_start, segment_end)
            print(segment.query_sequence)
            
            raise e
        
        if not seq_bases:  # Skip if no sequence was extracted
            continue
        
        metrics["segments_overlapping_region"] += 1
        
        # Convert to string
        sequence = ''.join(seq_bases)
        
        # Convert quality scores to FASTQ format (ASCII + 33)
        if qual_scores:
            quality = ''.join(chr(q + 33) for q in qual_scores)
        else:
            quality = '!' * len(sequence)  # Default to lowest quality if not available
        
        # Create read name with region information
        read_name = f"{segment.query_name}_{region_name}_{chrom}:{region_start}-{region_end}"
        
        # Write FASTQ entry
        fastq_out.write(f"@{read_name}\n")
        fastq_out.write(f"{sequence}\n")
        fastq_out.write(f"+\n")
        fastq_out.write(f"{quality}\n")

# Write metrics to YAML file if requested
if output_metrics:
    with open(output_metrics, 'w') as f:
        f.write("# Slice metrics\n")
        
        yaml.dump(metrics, f, default_flow_style=False) 