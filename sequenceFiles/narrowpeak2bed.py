#!/usr/bin/python
"""
Converts narrowPeak format to BED format.  

Usage: 

    narrowpeak2bed.py input.narrowbed > output.bed
"""
import sys
fn = sys.argv[1]
for line in open(fn):
    L = line.strip().split('\t')
    chrom = L[0]
    if not chrom.startswith('chr'):
        chrom = 'chr'+chrom
    newline = [chrom,L[1],L[2]]
    sys.stdout.write('\t'.join(newline)+'\n')

