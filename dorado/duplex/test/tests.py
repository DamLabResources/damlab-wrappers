import pytest
import os
import pysam

def test_output_exists():
    """Test that BAM file was created"""
    assert os.path.exists('test.bam')

def test_bam_is_valid():
    """Test that output is a valid BAM file"""
    try:
        bam = pysam.AlignmentFile('test.bam', 'rb', check_sq=False)
        # Try to read a record to verify it's a valid BAM
        next(bam)
        bam.close()
    except Exception as e:
        pytest.fail(f"BAM file is not valid: {str(e)}")

def test_bam_has_reads():
    """Test that BAM contains at least one read"""
    bam = pysam.AlignmentFile('test.bam', 'rb', check_sq=False)
    read_count = sum(1 for _ in bam)
    bam.close()
    assert read_count > 0, "BAM file contains no reads" 