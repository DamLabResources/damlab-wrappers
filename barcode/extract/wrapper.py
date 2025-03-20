"""Wrapper for extracting barcodes and UMIs from BAM files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"

import pysam
import regex
from contextlib import contextmanager

if "snakemake" not in locals():
    import snakemake # type: ignore

# Built-in regex patterns
BUILTIN_PATTERNS = {
    'SIVmac239m2': {
        'left_primer': '(CAAGCAGAAGACGGCATACGAGAT){i<=2,d<=2,s<=4}(?<umi>.{18})(atatacttagaaaaggaagaaggcatcataccaga){i<=2,d<=2,s<=4}',
        'right_primer': '(gaagaactccgaaaaaggctaaggc){i<=2,d<=2,s<=4}(?<umi>.{18})(GATCTCGGTGGTCGCCGTATCATT){i<=2,d<=2,s<=4}',
        'barcode_primer': '(cctcctccccctccaggactagcataa){i<=2,d<=2,s<=4}(?<umi>.{34})(atggaagaaagacctccagaaaatga){i<=2,d<=2,s<=4}'
    }
}

def get_regex_patterns(params):
    """Get regex patterns from either builtin or custom parameters"""
    if 'builtin' in params:
        if params['builtin'] not in BUILTIN_PATTERNS:
            raise ValueError(f"Unknown builtin pattern: {params['builtin']}. Available patterns: {list(BUILTIN_PATTERNS.keys())}")
        patterns = BUILTIN_PATTERNS[params['builtin']]
    else:
        required = ['left_primer', 'right_primer', 'barcode_primer']
        missing = [p for p in required if p not in params]
        if missing:
            raise ValueError(f"Missing required primer patterns: {missing}")
        patterns = {
            'left_primer': params['left_primer'],
            'right_primer': params['right_primer'],
            'barcode_primer': params['barcode_primer']
        }
    
    return {
        'LEFT_PRIMER': regex.compile(patterns['left_primer'], flags=regex.BESTMATCH+regex.IGNORECASE),
        'RIGHT_PRIMER': regex.compile(patterns['right_primer'], flags=regex.BESTMATCH+regex.IGNORECASE),
        'BARCODE_PRIMER': regex.compile(patterns['barcode_primer'], flags=regex.BESTMATCH+regex.IGNORECASE)
    }

def extract_umi(reg, text, default=None):
    """Extract UMI sequence from text using regex pattern"""
    matches = reg.finditer(text)
    for match in matches:
        return match.capturesdict()['umi'][0]
    return default

@contextmanager
def make_new_bamfile(insam, outpath):
    """Context manager for creating new BAM file"""
    with open(outpath, mode='wb') as handle:
        outsam = pysam.AlignmentFile(handle, mode='wb', header=insam.header)
        yield outsam

# Get input/output paths
bam_in_path = str(snakemake.input[0])
bam_out_path = str(snakemake.output[0])

# Get parameters with defaults
bc_tag = snakemake.params.get('barcode_tag', 'CR')  # Cellular barcode sequence
umi_tag = snakemake.params.get('umi_tag', 'OX')  # Original UMI bases

# Get regex patterns
patterns = get_regex_patterns(snakemake.params)

# Process BAM file
with open(bam_in_path, mode='rb') as handle:
    samfile = pysam.AlignmentFile(handle)
    with make_new_bamfile(samfile, bam_out_path) as newsam:
        for read in samfile:
            if read.query_sequence:
                # Extract barcode
                barcode = extract_umi(patterns['BARCODE_PRIMER'], read.query_sequence, default='N'*34)
                read.set_tag(bc_tag, barcode, 'Z')
                
                # Extract UMIs
                left_umi = extract_umi(patterns['LEFT_PRIMER'], read.query_sequence, default='N'*18)
                right_umi = extract_umi(patterns['RIGHT_PRIMER'], read.query_sequence, default='N'*18)
                read.set_tag(umi_tag, left_umi+right_umi, 'Z')
                
                newsam.write(read)
