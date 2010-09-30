#!/usr/bin/python
"""
Reads stdin, writes stdout.  Converts a SAM stream into a naive Bowtie outpuf
format (that is, it ignores CIGAR strings)

Typical usage:

    samtools view -S -F 0x0004 $SAMFILE | python sam2bowtie.py > treatment.sam
"""
import sys
for line in sys.stdin:
    if line.startswith('@'):
        continue
    L = line.split()
    strand = L[1]
    if strand == '0':
        strand = '+'
    elif strand == '16':
        strand = '-'
    else:
        continue
    chrom = L[2]
    start = L[3]
    name = L[0]
    seq = L[9]
    qual = L[10]
    new = [name,strand,chrom,start,seq,qual,'0','']
    sys.stdout.write('\t'.join(new)+'\n')
