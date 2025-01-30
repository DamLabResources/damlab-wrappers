import pytest
import os
from Bio import SeqIO  # type: ignore

def test_haplotypes_output_exists():
    """Test that haplotypes FASTA file was created"""
    assert os.path.exists('scratch/test_haplotypes.fasta')

def test_directory_output_exists():
    """Test that output directory was created"""
    assert os.path.exists('scratch/test_directory')
    assert os.path.exists('scratch/test_directory/haplotypes.final.fa')

def test_haplotypes_is_valid_fasta():
    """Test that haplotypes output is a valid FASTA file"""
    try:
        with open('scratch/test_haplotypes.fasta', 'r') as f:
            records = list(SeqIO.parse(f, 'fasta'))
        assert len(records) > 0, "FASTA file contains no sequences"
    except Exception as e:
        pytest.fail(f"FASTA file is not valid: {str(e)}")

def test_directory_haplotypes_is_valid_fasta():
    """Test that directory output contains valid FASTA file"""
    try:
        with open('scratch/test_directory/haplotypes.final.fa', 'r') as f:
            records = list(SeqIO.parse(f, 'fasta'))
        assert len(records) > 0, "Directory FASTA file contains no sequences"
    except Exception as e:
        pytest.fail(f"Directory FASTA file is not valid: {str(e)}")

def test_haplotypes_have_reasonable_length():
    """Test that haplotype sequences are of reasonable length"""
    with open('scratch/test_haplotypes.fasta', 'r') as f:
        records = list(SeqIO.parse(f, 'fasta'))
        for record in records:
            assert len(record.seq) > 50, "Haplotype sequence is too short"
            assert len(record.seq) < 100000, "Haplotype sequence is too long"

def test_directory_contains_expected_files():
    """Test that output directory contains expected auxiliary files"""
    expected_files = [
        'haplotypes.final.fa',
        'reads.fasta'
    ]
    for file in expected_files:
        assert os.path.exists(os.path.join('scratch/test_directory', file)), \
            f"Expected file {file} missing from output directory" 