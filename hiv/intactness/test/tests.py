import os
import pytest
import yaml

def test_outputs_exist():
    """Test that output YAML files were created"""
    assert os.path.exists('test_fasta_output.yaml')
    assert os.path.exists('test_sam_output.yaml')

def check_yaml_format(stats):
    """Helper function to check YAML format"""
    required_keys = {
        'total_reads',
        'countable_reads',
        'intact_reads',
        'percent_countable',
        'percent_intact',
        'sample_name'
    }
    
    assert set(stats.keys()) == required_keys
    assert isinstance(stats['total_reads'], int)
    assert isinstance(stats['countable_reads'], int)
    assert isinstance(stats['intact_reads'], int)
    assert isinstance(stats['percent_countable'], (int, float))
    assert isinstance(stats['percent_intact'], (int, float))
    assert isinstance(stats['sample_name'], str)

def check_counts(stats):
    """Helper function to check count logic"""
    assert stats['countable_reads'] <= stats['total_reads']
    assert stats['intact_reads'] <= stats['countable_reads']
    assert 0 <= stats['percent_countable'] <= 100
    assert 0 <= stats['percent_intact'] <= 100

def test_fasta_output():
    """Test FASTA output format and counts"""
    with open('test_fasta_output.yaml') as f:
        stats = yaml.safe_load(f)
    check_yaml_format(stats)
    check_counts(stats)
    
    # Test specific counts for FASTA
    assert stats['total_reads'] == 6
    assert stats['countable_reads'] == 4  # 60, 66, 72, 78 bp sequences
    assert stats['intact_reads'] == 1     # 72 bp sequences

def test_sam_output():
    """Test SAM output format and counts"""
    with open('test_sam_output.yaml') as f:
        stats = yaml.safe_load(f)
    check_yaml_format(stats)
    check_counts(stats)
    
    # Test specific counts for SAM
    assert stats['total_reads'] == 6
    assert stats['countable_reads'] == 4  # 60M, 66M, 72M, 78M alignments
    assert stats['intact_reads'] == 1     # 72M alignments 