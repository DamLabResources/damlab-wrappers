"""Unit tests for Biopython translate wrapper"""

import os
import pytest
from Bio import SeqIO

def test_output_exists():
    """Test that output files were created"""
    assert os.path.exists('test_output/basic.fasta'), "Basic output not found"
    assert os.path.exists('test_output/frame1.fasta'), "Frame 1 output not found"
    
    assert os.path.exists('test_output/basic.log'), "Basic log not found"
    assert os.path.exists('test_output/frame1.log'), "Frame 1 log not found"

def test_log_content():
    """Test that log files contain expected content"""
    with open('test_output/basic.log', 'r') as log_file:
        assert "Processed 2 sequences" in log_file.read(), "Basic log does not contain expected content"

    with open('test_output/frame1.log', 'r') as log_file:
        assert "Processed 2 sequences" in log_file.read(), "Frame 1 log does not contain expected content"


def test_basic_translation():
    """Test basic translation with default parameters"""
    records = list(SeqIO.parse('test_output/basic.fasta', 'fasta'))
    assert len(records) == 2, "Expected 2 sequences in output"
    
    # Check first sequence
    assert records[0].id == "seq1", "Expected sequence ID to be preserved"
    assert str(records[0].seq) == "MAMAK", "Incorrect translation"
    
    # Check second sequence
    assert records[1].id == "seq2", "Expected sequence ID to be preserved"
    assert str(records[1].seq) == "MAMAK", "Incorrect translation"

def test_frame1_translation():
    """Test translation with frame 1"""
    records = list(SeqIO.parse('test_output/frame1.fasta', 'fasta'))
    assert len(records) == 2, "Expected 2 sequences in output"
    
    # Check first sequence
    assert records[0].id == "seq1", "Expected sequence ID to be preserved"
    assert str(records[0].seq) == "WPWP", "Incorrect translation in frame 1"
    
    # Check second sequence
    assert records[1].id == "seq2", "Expected sequence ID to be preserved"
    assert str(records[1].seq) == "WPWP", "Incorrect translation in frame 1"

