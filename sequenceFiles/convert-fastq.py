#!/usr/bin/python

usage = """
Script to convert FASTQ files between quality score formats. Requires BioPython v1.5+

Usage:

    python convert-to-illumina.py --informat fastq-sanger --outformat fastq-illumina SRX00001.fastq > SRX00001.illumina.fastq

--informat and --outformat are those formats supported by BioPython's SeqIO. Useful ones include:

fastq-sanger
fastq-solexa
fastq-illumina
"""
import optparse
import sys
from Bio import SeqIO

op = optparse.OptionParser(usage=usage)
op.add_option('--informat',dest='informat',help='Input format')
op.add_option('--outformat',dest='outformat',help='Output format')
options,args = op.parse_args()
if len(args) != 1:
    op.print_help()
    print 'Need exactly one input file'
    sys.exit(1)

if (not options.informat) or (not options.outformat):
    op.print_help()
    print 'Need to specify both --informat and --outformat'
    sys.exit(1)
    
fn = args[0]
for rec in SeqIO.parse(open(fn),options.informat):
    sys.stdout.write(rec.format(options.outformat))
    
