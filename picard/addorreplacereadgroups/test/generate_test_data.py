"""Generate test BAM file for AddOrReplaceReadGroups wrapper tests"""

import pysam
import os

def create_test_bam(filename):
    """Create a simple BAM file with a few reads"""
    
    # Make sure output directory exists
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    
    # Create header
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 1000, 'SN': 'chr1'}]
    }
    
    # Create BAM file
    with pysam.AlignmentFile(filename, 'wb', header=header) as outf:
        # Create a few dummy reads
        for i in range(10):
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i}"
            a.query_sequence = "ACGT" * 10
            a.flag = 0
            a.reference_id = 0
            a.reference_start = i * 100
            a.mapping_quality = 20
            a.cigar = ((0, 40),)
            a.query_qualities = [30] * 40
            outf.write(a)

if __name__ == "__main__":
    create_test_bam("test_data/input.bam") 