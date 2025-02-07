import pytest # type: ignore
from pathlib import Path
import subprocess

def test_pod5_output_exists():
    """Test POD5 file exists"""
    assert Path('test_single.pod5').exists()
    assert Path('test_dir.pod5').exists()


def check_pod5_file_is_valid(pod5_file: Path):
    """Check POD5 file is valid"""
    result = subprocess.run(['pod5', 'view', pod5_file], 
                          capture_output=True, 
                          text=True)
        
    assert result.returncode == 0
    assert len(result.stdout.split('\n')) == 6, 'Expected 6 lines of output from POD5 file'

def test_pod5_is_valid():
    """Test POD5 file is valid"""
    check_pod5_file_is_valid(Path('test_single.pod5'))

def test_pod5_dir_is_valid():
    """Test POD5 directory is valid"""
    check_pod5_file_is_valid(Path('test_dir.pod5'))

