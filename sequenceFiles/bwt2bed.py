#!/usr/bin/python

import sys

usage = """Converts a Bowtie output file to a BED file.  BED feature names are
taken from the first column. Scores are set to zero.

Example usage:

    python bwt2bed.py myreads.bwt > myreads.bed"""
try:
    bwt = sys.argv[1]
except IndexError:
    print usage
    sys.exit(1)


for line in open(bwt):
    L = line.split()
    name = L[0]
    strand = L[1]
    chrom = L[2]
    start = L[3] 
    seq = L[4]
    score = '0'
    stop = str(int(start)+len(seq))

    bedline = [chrom,start,stop,name,score,strand]
    sys.stdout.write('\t'.join(bedline)+'\n')


