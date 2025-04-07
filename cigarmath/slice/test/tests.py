"""Unit tests for slice functions"""

import os
import pytest
import yaml

def test_output_exists():
    """Test that output files were created"""
    assert os.path.exists('test_output.fastq')
    assert os.path.exists('test_metrics.yaml')

def test_output_fastq_format():
    """Test that output FASTQ has correct format"""
    with open('test_output.fastq') as f:
        lines = f.readlines()
        
    # Check if file is not empty
    assert len(lines) > 0
    
    # Check FASTQ format (groups of 4 lines)
    assert len(lines) % 4 == 0
    
    for i in range(0, len(lines), 4):
        # Line 1: Header starts with @
        assert lines[i].startswith('@')
        # Line 3: Plus line
        assert lines[i+2].strip() == '+'
        # Line 4: Quality score length matches sequence length
        assert len(lines[i+1].strip()) == len(lines[i+3].strip())
        
        # Check that read name contains region info
        assert 'HXB2:110-130' in lines[i] or 'test_region' in lines[i]

def test_metrics_yaml():
    """Test that metrics YAML contains expected fields"""
    with open('test_metrics.yaml') as f:
        metrics = yaml.safe_load(f)
    
    # Check required fields
    assert 'region' in metrics
    assert 'region_name' in metrics
    assert 'sample_name' in metrics
    assert 'total_segments_processed' in metrics
    assert 'segments_overlapping_region' in metrics
    
    # Verify values
    assert metrics['region'] == 'HXB2:110-130'
    assert metrics['region_name'] == 'test_region'
    assert metrics['sample_name'] == 'test_sample'
    assert metrics['total_segments_processed'] > 0
    assert metrics['segments_overlapping_region'] > 0
    
    # Verify that we found at least 3 overlapping segments
    # (should be reads 3, 4, and 5 from our test SAM)
    assert metrics['segments_overlapping_region'] >= 3

def test_read_content():
    """Test that extracted reads contain expected content"""
    with open('test_output.fastq') as f:
        content = f.read()
    
    # Check for presence of specific sequences from overlapping reads
    assert 'GTACGTACGT' in content  # Part of read3
    assert 'TACGTAC' in content     # Part of read4
    
    # Make sure we don't have reads outside the region
    with open('test_metrics.yaml') as f:
        metrics = yaml.safe_load(f)
    
    # We requested region 110-130, so read6 (position 200) shouldn't be included
    assert metrics['total_segments_processed'] > metrics['segments_overlapping_region'] 