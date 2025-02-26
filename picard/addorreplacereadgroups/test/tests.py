import os
import pysam # type: ignore
import pytest # type: ignore

def test_output_files_exist():
    """Test that output BAM files were created"""
    assert os.path.exists('test_output/basic.bam'), "Basic output BAM not found"
    assert os.path.exists('test_output/full.bam'), "Full output BAM not found"
    assert os.path.exists('test_output/minimal.bam'), "Minimal output BAM not found"

def test_output_files_valid():
    """Test that output BAM files are valid"""
    for bam_file in ['test_output/basic.bam', 'test_output/full.bam', 'test_output/minimal.bam']:
        try:
            bam = pysam.AlignmentFile(bam_file, 'rb')
            next(bam)
            bam.close()
        except Exception as e:
            pytest.fail(f"BAM file {bam_file} is not valid: {str(e)}")

def test_basic_read_groups():
    """Test basic read group parameters are set correctly"""
    bam = pysam.AlignmentFile('test_output/basic.bam', 'rb')
    header = bam.header
    
    assert 'RG' in header, "No read groups in header"
    rg = header['RG'][0]  # Get first read group
    
    assert rg['ID'] == 'test_id'
    assert rg['LB'] == 'test_lib'
    assert rg['PL'] == 'ILLUMINA'
    assert rg['PU'] == 'unit1'
    assert rg['SM'] == 'sample1'
    
    bam.close()

def test_full_read_groups():
    """Test all read group parameters are set correctly"""
    bam = pysam.AlignmentFile('test_output/full.bam', 'rb')
    header = bam.header
    
    assert 'RG' in header, "No read groups in header"
    rg = header['RG'][0]  # Get first read group
    
    # Check required parameters
    assert rg['ID'] == 'full_id'
    assert rg['LB'] == 'full_lib'
    assert rg['PL'] == 'ILLUMINA'
    assert rg['PU'] == 'unit2'
    assert rg['SM'] == 'sample2'
    
    # Check optional parameters
    assert rg['CN'] == 'test_center'
    assert rg['DS'] == 'test description'
    assert rg['PM'] == 'NextSeq2000'
    
    bam.close()

def test_minimal_read_groups():
    """Test minimal read group parameters with ONT platform"""
    bam = pysam.AlignmentFile('test_output/minimal.bam', 'rb')
    header = bam.header
    
    assert 'RG' in header, "No read groups in header"
    rg = header['RG'][0]  # Get first read group
    
    assert rg['ID'] == 'min_id'
    assert rg['LB'] == 'min_lib'
    assert rg['PL'] == 'ONT'
    assert rg['PU'] == 'unit3'
    assert rg['SM'] == 'sample3'
    
    bam.close()

def test_read_groups_assigned():
    """Test that reads in output BAMs have read groups assigned"""
    for bam_file in ['test_output/basic.bam', 'test_output/full.bam', 'test_output/minimal.bam']:
        bam = pysam.AlignmentFile(bam_file, 'rb')
        read = next(bam)
        assert read.get_tag('RG'), f"Read in {bam_file} missing RG tag"
        bam.close() 