import pytest
import os
import yaml
import csv

def test_output_exists():
    """Test that output files exist"""
    assert os.path.exists('results.csv')
    assert os.path.exists('summary.yaml')

def test_csv_format():
    """Test CSV file format and content"""
    with open('results.csv', 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        
        # Check header format
        assert header[0] == 'read_id'
        assert sorted(header[1:]) == sorted(['p1', 'p2', 'p3', 'p4', 'P5', 'p6'])
        
        # Convert rows to dict for easier testing
        results = {row[0]: {h: v for h, v in zip(header[1:], row[1:])} 
                  for row in reader}
        
        # Check specific known results
        assert results['seq1']['p1'] != 'None'  # Should have a length
        assert results['seq1']['p2'] != 'None'  # Should have a length
        assert results['seq1']['p3'] != 'None'  # Should have a length
        assert results['seq2']['p4'] != 'None'  # Should have a length
        assert results['seq2']['P5'] != 'None'  # Should have a length
        assert results['seq2']['p6'] != 'None'  # Should have a length

def test_summary_format():
    """Test YAML summary format and content"""
    with open('summary.yaml', 'r') as f:
        summary = yaml.safe_load(f)
        
        # Check structure
        assert 'total_sequences' in summary
        assert 'primer_hits' in summary
        assert 'pairwise_hits' in summary
        
        # Check content
        assert summary['total_sequences'] == 2  # seq1 and seq2
        
        # Check primer hits
        assert summary['primer_hits']['p1'] == 1  # Only hits seq1
        assert summary['primer_hits']['p4'] == 2  # Hits both sequences

        # Check sample name
        assert summary['sample_name'] == 'sample1'
        