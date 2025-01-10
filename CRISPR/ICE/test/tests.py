import pytest
import os
import json

# Define tests here that validate rules
# Add test-specific libraries to env.yaml

# At a minimum, tests should validate that test rules have created
# files in the relevant formats, applicable headers, etc.


def check_paths_equal(path1, path2):
    with open(path1) as h1:
        with open(path2) as h2:
            if 'json' in path1:
                assert json.load(h1) == json.load(h2)
            else:
                assert h1.read() == h2.read()

def test_all_json():
    check_paths_equal('check.all.json', 'testing.all.json')
    
def test_all_txt():
    check_paths_equal('check.all.txt', 'testing.all.txt')

def test_contribs_json():
    check_paths_equal('check.contribs.json', 'testing.contribs.json')

def test_contribs_txt():
    check_paths_equal('check.contribs.txt', 'testing.contribs.txt')

def test_indel_json():
    check_paths_equal('check.indel.json', 'testing.indel.json')
    
def test_trace_json():
    check_paths_equal('check.trace.json', 'testing.trace.json')

def test_windowed_json():
    check_paths_equal('check.windowed.json', 'testing.windowed.json')
    
def test_windowed_txt():
    check_paths_equal('check.windowed.txt', 'testing.windowed.txt')
