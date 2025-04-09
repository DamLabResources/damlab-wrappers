"""Wrapper for correcting barcodes in BAM files using UMI-tools"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.0"

import pysam
from contextlib import contextmanager
from collections import Counter
from umi_tools import UMIClusterer
import yaml
from typing import Dict, Any

if "snakemake" not in locals():
    import snakemake # type: ignore

@contextmanager
def make_new_bamfile(insam, outpath):
    """Context manager for creating new BAM file"""
    with open(outpath, mode='wb') as handle:
        header = insam.header.to_dict()
        
        # Add program group to header
        if 'PG' not in header:
            header['PG'] = []
            
        header['PG'].append({
            'ID': 'barcode_correct',
            'PN': 'damlab-barcode-correct',
            'VN': __version__,
            'DS': 'Correct barcodes using UMI-tools clustering'
        })
        
        outsam = pysam.AlignmentFile(handle, mode='wb', header=header)
        yield outsam

def get_counts(samfile, tag, bc_length=18):
    """Takes a samfile and tag and returns a count of each barcode"""
    counter = Counter()
    for num, read in enumerate(samfile):
        try:
            barcode = read.get_tag(tag)
        except KeyError:
            continue
        counter[bytes(barcode, 'ascii')] += 1
    return counter

def make_clusters(counter_dict, mismatches=3):
    """With a {umi:count} dict, return a {original:corrected} dict"""
    clusterer = UMIClusterer(cluster_method="adjacency")
    clusters = clusterer(counter_dict, mismatches)
    
    cluster_dict = {}
    for cluster in clusters:
        key = cluster[0]
        for item in cluster:
            cluster_dict[item] = key
    
    return cluster_dict

def correct_tag_for_stream(stream, in_tag, out_tag, cluster_mapping, bc_length=18):
    """Yield reads with corrected tags"""
    for read in stream:
        try:
            bc = read.get_tag(in_tag)
        except KeyError:
            continue
        if len(bc) == bc_length:
            corrected = cluster_mapping[bytes(bc,'ascii')]
            read.set_tag(out_tag, corrected, 'Z')
            yield read

def generate_metrics(
    original_counts: Dict[bytes, int],
    cluster_mapping: Dict[bytes, bytes],
    corrected_counts: Dict[bytes, int]
) -> Dict[str, Any]:
    """Generate metrics from barcode correction results.
    
    Args:
        original_counts: Dictionary of original barcode counts
        cluster_mapping: Dictionary mapping original to corrected barcodes
        corrected_counts: Dictionary of counts after correction
        
    Returns:
        Dictionary of metrics including counts and histograms
    """
    # Convert bytes keys to strings for YAML serialization
    orig_counts_str = {k.decode('ascii'): v for k, v in original_counts.items()}
    corr_counts_str = {k.decode('ascii'): v for k, v in corrected_counts.items()}
    
    # Generate count histogram for corrected barcodes
    count_hist = {}
    for count in corrected_counts.values():
        count_hist[count] = count_hist.get(count, 0) + 1
    
    metrics = {
        'barcode_counts': {
            'unique_barcodes_before': len(original_counts),
            'unique_barcodes_after': len(corrected_counts),
            'total_reads': sum(original_counts.values())
        },
        'correction_details': {
            'clusters_formed': len(set(cluster_mapping.values())),
            'barcodes_corrected': len(cluster_mapping) - len(corrected_counts)
        },
        'count_data': {
            'original_counts': orig_counts_str,
            'corrected_counts': corr_counts_str,
            'count_histogram': dict(sorted(count_hist.items()))
        }
    }
    
    return metrics

# Get input/output paths
bam_in_path = str(snakemake.input[0])
bam_out_path = str(snakemake.output[0])

# Get parameters
in_tag = snakemake.params['in_tag']
out_tag = snakemake.params['out_tag']
bc_len = snakemake.params['barcode_length']
mismatches = snakemake.params.get('mismatches', 3)
sample_name = snakemake.params.get('sample_name', None)
barcode_name = snakemake.params.get('barcode_name', None)

# Process BAM file
with open(bam_in_path, mode='rb') as handle:
    insam = pysam.AlignmentFile(handle, mode='rb')
    barcode_counts = get_counts(insam, in_tag, bc_length=bc_len)

cluster_mapping = make_clusters(barcode_counts, mismatches=mismatches)

# Calculate corrected counts
corrected_counts = {}
for orig_bc, count in barcode_counts.items():
    corrected = cluster_mapping[orig_bc]
    corrected_counts[corrected] = corrected_counts.get(corrected, 0) + count

# Get metrics file path if specified
metrics_file = snakemake.output.get('metrics', None)

with open(bam_in_path, mode='rb') as handle:
    insam = pysam.AlignmentFile(handle, mode='rb')
    with make_new_bamfile(insam, bam_out_path) as newsam:
        stream = correct_tag_for_stream(insam, in_tag, out_tag, 
                                      cluster_mapping,
                                      bc_length=bc_len)
        for read in stream:
            newsam.write(read)

# Write metrics if requested
if metrics_file:
    metrics = generate_metrics(barcode_counts, cluster_mapping, corrected_counts)
    if sample_name:
        metrics['sample_name'] = sample_name
    if barcode_name:
        metrics['barcode_name'] = barcode_name
    with open(metrics_file, 'w') as f:
        f.write('# Barcode correction metrics\n')
        yaml.dump(metrics, f, default_flow_style=False) 