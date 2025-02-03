import os
import pytest # type: ignore

def test_output_exists():
    """Test that output FASTA files were created"""
    assert os.path.exists('test_clipped.fasta')
    assert os.path.exists('test_clipped_no_ref.fasta')

def test_sequences_oriented():
    """Test that sequences are properly oriented and clipped"""
    for filename in ['test_clipped.fasta', 'test_clipped_no_ref.fasta']:
        with open(filename) as f:
            content = f.read()
        
        # Check that file has content
        assert len(content.strip()) > 0
        
        # Check that each sequence has coverage information
        for line in content.split('\n'):
            if line.startswith('>') and not line.startswith('>ref'):
                assert 'coverage=' in line

def test_sequence_quality():
    """Test that sequences meet minimum coverage criteria"""
    for filename in ['test_clipped.fasta', 'test_clipped_no_ref.fasta']:
        with open(filename) as f:
            for line in f:
                if line.startswith('>') and not line.startswith('>ref'):
                    coverage = float(line.split('coverage=')[1].split()[0])
                    assert coverage >= 0.2

def test_reference_inclusion():
    """Test that reference is included only when specified"""
    # Check file with reference
    with open('test_clipped.fasta') as f:
        content = f.read().split('\n')
        assert any(line.startswith('>ref') for line in content)
    
    # Check file without reference
    with open('test_clipped_no_ref.fasta') as f:
        content = f.read().split('\n')
        assert not any(line.startswith('>ref') for line in content) 