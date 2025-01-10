__author__ = "Will Dampier"
__copyright__ = "Copyright 2022, Will Dampier"
__email__ = "wnd22@drexel.edu"
__license__ = "Kept"

"""Snakemake wrapper for synthego_ice."""

from os import path
from shutil import copy
from snakemake.shell import shell
from tempfile import TemporaryDirectory


extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

assert snakemake.params.get('target'), "Must provide a target parameter!"

with TemporaryDirectory(delete=False) as path:
    
    cmd = [f"docker run -it -v {path}:/data -w /ice",
           '--user "$(id -u):$(id -g)"',
           "-i synthego/ice",
           "python ice_analysis_single.py",
           "--control /data/control.ab1",
           "--edited /data/edited.ab1",
           "--target {snakemake.params.target}",
           "--out /data/results/testing"]
    
    if 'donor' in snakemake.params:
        cmd.append("--donor {snakemake.params.donor}")
    
    copy(snakemake.input.control, path+'/control.ab1')
    copy(snakemake.input.edited, path+'/edited.ab1')
        
    shell(' '.join(cmd))
    
    exts = ['all.json',
            'all.txt',
            'contribs.json',
            'contribs.txt',
            'indel.json',
            'trace.json',
            'windowed.json',
            'windowed.txt']
    
    for ext in exts:
        key = ext.replace('.', '_')
        #print('checking for', key)
        if snakemake.output.get(key):
            #print('found')
            copy(path + '/results/testing.'+ext, snakemake.output[key])