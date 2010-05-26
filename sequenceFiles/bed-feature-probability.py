#!/usr/bin/python

import optparse
import sys
import logging

usage = """

\t%prog -i BEDFILE [--genomesize=GENOMESIZE]

Script to calculate the probability of encountering a feature in the input bed
file assuming random distribution and random sampling.  Output is to stdout.
"""
op = optparse.OptionParser(usage=usage)
op.add_option('-i',dest='bedfile',
              help='Input BED format file. Only uses first 3 fields.')
op.add_option('--genomesize', dest='genomesize',default=120e6,
              help='Size of genome to compare against; default=%default')
options,args = op.parse_args()
if options.bedfile is None:
    logging.warning('Need an input file. Use -h for help.')
    sys.exit()

total_coverage = 0
for line in open(options.bedfile):
    if 'track' in line or 'browser' in line:
        continue
    L = line.rstrip().split('\t')
    start = int(L[1])
    stop = int(L[2])
    length = abs(start-stop)
    total_coverage += length

prob = total_coverage / float(options.genomesize)
print """

Total coverage of sequences in %s:
    
    %s

Probability of randomly encountering a sequence in %s
in a genome of size %d: 
    
    %.2f""" % (options.bedfile,total_coverage, options.bedfile,options.genomesize,prob)
