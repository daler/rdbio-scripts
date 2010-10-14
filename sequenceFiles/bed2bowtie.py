#!/usr/bin/python
"""
Converts a BED file into a Bowtie format file.  Since Bowtie files contain much
more information than a bed file, this script creates fake sequences and
quality scores.

Yep, it's a hack but some programs (e.g., spp) expect Bowtie format files.

Usage:

    bed2bowtie.py input.bed > output.bowtie
"""

import sys
fn = sys.argv[1]

for line in open(fn):
    L = line.strip().split('\t')
    chrom,start,stop,name,score,strand = L[0:6]
    seq = 'N'*(int(stop)-int(start))
    qual = '^'*len(seq)
    newline = ['FAKE',strand,chrom,start,seq,qual,'0','']
    sys.stdout.write('\t'.join(newline)+'\n')



