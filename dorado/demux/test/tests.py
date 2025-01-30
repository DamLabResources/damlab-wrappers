import pytest
import os
import pysam # type: ignore
from pathlib import Path

def test_output_dir_exists():
    """Test that output directory was created"""
    assert os.path.exists('test_bam')
    assert os.path.exists('test_fastq')

def test_bam_outputs_exist():
    """Test that BAM output files were created"""
    # Check for BAM outputs in single file test
    bam_files = list(Path('test_bam').glob('*.bam'))
    assert len(bam_files) > 0, "No BAM files found in single file output"


def test_fastq_outputs_exist():
    """Test that FASTQ output files were created"""
    # Check for FASTQ outputs in fastq test
    fastq_files = list(Path('test_fastq').glob('*.fastq'))
    assert len(fastq_files) > 0, "No FASTQ files found in output"

def test_bam_files_valid():
    """Test that output BAM files are valid"""
    for bam_file in Path('test_bam').glob('*.bam'):
        try:
            bam = pysam.AlignmentFile(str(bam_file), 'rb', check_sq=False)
            next(bam)
            bam.close()
        except Exception as e:
            pytest.fail(f"BAM file {bam_file} is not valid: {str(e)}")
            

def test_fastq_files_have_content():
    """Test that FASTQ files contain reads"""
    for fastq_file in Path('test_fastq').glob('*.fastq'):
        with open(fastq_file) as f:
            content = f.read()
            assert len(content) > 0, f"FASTQ file {fastq_file} is empty" 