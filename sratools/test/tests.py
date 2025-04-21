import os
import pytest
from Bio import SeqIO

def test_fastq_file_exists():
    """Test that the FASTQ file exists"""
    assert os.path.exists("test_output/DRR537798.fastq")

def test_fastq_format():
    """Test that the file is valid FASTQ format"""
    records = list(SeqIO.parse("test_output/DRR537798.fastq", "fastq"))
    assert len(records) > 10, "FASTQ file should have more than 10 records"
    
    # Check that these are single-read sequences
    for record in records:
        # Check that the ID doesn't end with /1 or /2 (which would indicate paired reads)
        assert not record.id.endswith("/1") and not record.id.endswith("/2"), \
            "Sequences should be single-read (no /1 or /2 suffixes)"

def test_fastq_quality():
    """Test that the FASTQ file has quality scores"""
    records = list(SeqIO.parse("test_output/DRR537798.fastq", "fastq"))
    for record in records:
        assert record.letter_annotations, "Records should have quality scores"
        assert len(record.letter_annotations["phred_quality"]) == len(record.seq), "Quality scores should match sequence length"
        
        # Check that quality scores are valid (between 0 and 93 for Illumina)
        qualities = record.letter_annotations["phred_quality"]
        assert all(0 <= q <= 93 for q in qualities), "Quality scores should be valid Illumina scores"
