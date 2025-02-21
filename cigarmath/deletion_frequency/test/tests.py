import os
import pytest # type: ignore
import yaml

def test_output_exists():
    """Test that output YAML file was created"""
    assert os.path.exists('test_output.yaml')

def test_output_format():
    """Test that output YAML has correct structure"""
    with open('test_output.yaml') as f:
        stats = yaml.safe_load(f)
    
    required_keys = {
        'sample_name',
        'total_reads',
        'reads_covering_required',
        'reads_with_deletion',
        'deletion_frequency'
    }
    
    assert set(stats.keys()) == required_keys
    
    correct_values = {
        'deletion_frequency': 13/30,
        'reads_covering_required': 30,
        'reads_with_deletion': 13,
        'sample_name': 'test',
        'total_reads': 180
    }

    for key, value in correct_values.items():
        if isinstance(value, str    ):
            assert stats[key] == value
        else:
            assert abs(stats[key] - value) < 0.001

