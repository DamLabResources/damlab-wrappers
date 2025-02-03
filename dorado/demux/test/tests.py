import os
from pathlib import Path
import pytest # type: ignore
import pysam # type: ignore


def test_bam_outputs_exist():
    """Test that BAM output files were created"""
    # Check for mapped output
    assert os.path.exists('test_bam/sample3.bam'), "Mapped BAM file not found"
    assert os.path.exists('test_bam/sample5.bam'), "Mapped BAM file not found"
    assert os.path.exists('test_bam/missing.bam'), "Empty BAM file not created"

def test_fastq_outputs_exist():
    """Test that FASTQ output files were created"""
    # Check for mapped output
    assert os.path.exists('test_fastq/sample3.fastq'), "Mapped FASTQ file not found"
    assert os.path.exists('test_fastq/sample5.fastq'), "Mapped FASTQ file not found"
    assert os.path.exists('test_fastq/missing.fastq'), "Empty FASTQ file not created"

def test_bam_files_valid():
    """Test that output BAM files are valid"""
    for bam_file in Path('test_bam').glob('**/*.bam'):
        if bam_file.name == 'missing.bam':
            # Check that missing.bam is empty
            assert os.path.getsize(bam_file) == 0, "Missing BAM file should be empty"
            continue
            
        try:
            bam = pysam.AlignmentFile(str(bam_file), 'rb', check_sq=False)
            next(bam)
            bam.close()
        except Exception as e:
            pytest.fail(f"BAM file {bam_file} is not valid: {str(e)}")

def test_fastq_files_have_content():
    """Test that FASTQ files contain reads"""
    for fastq_file in Path('test_fastq').glob('**/*.fastq'):
        if fastq_file.name == 'missing.fastq':
            # Check that missing.fastq is empty
            assert os.path.getsize(fastq_file) == 0, "Missing FASTQ file should be empty"
            continue
            
        with open(fastq_file) as f:
            content = f.read()
            assert len(content) > 0, f"FASTQ file {fastq_file} is empty" 