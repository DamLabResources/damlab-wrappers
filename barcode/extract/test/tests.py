import pysam
import os

def test_output_exists():
    """Test that output BAM file was created"""
    assert os.path.exists('test_output/test_output.bam')

def test_output_tags():
    """Test that output BAM has correct tags"""
    with pysam.AlignmentFile('test_output/test_output.bam', 'rb') as bam:
        for read in bam:
            if read.query_sequence:
                assert read.has_tag('CR')  # Has barcode tag
                assert read.has_tag('OX')  # Has UMI tag
                
                # Check tag lengths
                assert len(read.get_tag('CR')) == 34  # Barcode length
                assert len(read.get_tag('OX')) == 36  # Combined UMI length 