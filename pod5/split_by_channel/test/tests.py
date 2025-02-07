import pytest # type: ignore
from pathlib import Path
import subprocess
import csv

def test_output_dirs_exist():
    """Test output directories exist"""
    assert Path('test_single_split').exists()
    assert Path('test_dir_split').exists()


def check_pod5_dir_is_valid(pod5_dir: Path):
    """Check POD5 directory contains valid files"""
    # Should have at least one POD5 file
    pod5_files = list(pod5_dir.glob("*.pod5"))
    assert len(pod5_files) > 0, f"No POD5 files found in {pod5_dir}"

    # Check each file is valid
    for pod5_file in pod5_files:
        result = subprocess.run(['pod5', 'view', pod5_file],
                              capture_output=True,
                              text=True)
        assert result.returncode == 0
        assert len(result.stdout.split('\n')) > 1, f'No reads found in {pod5_file}'

def test_split_outputs_are_valid():
    """Test split outputs are valid"""
    check_pod5_dir_is_valid(Path('test_single_split'))
    check_pod5_dir_is_valid(Path('test_dir_split'))
