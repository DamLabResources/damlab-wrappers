"""Wrapper for correcting barcodes in BAM files using UMI-tools"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import pysam
from contextlib import contextmanager
from collections import Counter
from umi_tools import UMIClusterer

if "snakemake" not in locals():
    import snakemake # type: ignore

@contextmanager
def make_new_bamfile(insam, outpath):
    """Context manager for creating new BAM file"""
    with open(outpath, mode='wb') as handle:
        outsam = pysam.AlignmentFile(handle, mode='wb', header=insam.header)
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

# Get input/output paths
bam_in_path = str(snakemake.input[0])
bam_out_path = str(snakemake.output[0])

# Get parameters
in_tag = snakemake.params['in_tag']
out_tag = snakemake.params['out_tag']
bc_len = snakemake.params['barcode_length']
mismatches = snakemake.params.get('mismatches', 3)

# Process BAM file
with open(bam_in_path, mode='rb') as handle:
    insam = pysam.AlignmentFile(handle, mode='rb')
    barcode_counts = get_counts(insam, in_tag, bc_length=bc_len)

cluster_mapping = make_clusters(barcode_counts, mismatches=mismatches)

with open(bam_in_path, mode='rb') as handle:
    insam = pysam.AlignmentFile(handle, mode='rb')
    with make_new_bamfile(insam, bam_out_path) as newsam:
        stream = correct_tag_for_stream(insam, in_tag, out_tag, 
                                      cluster_mapping,
                                      bc_length=bc_len)
        for read in stream:
            newsam.write(read) 