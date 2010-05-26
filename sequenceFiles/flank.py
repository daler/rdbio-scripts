#!/usr/bin/python

__doc__ = """
Script to print the flanking regions of a bed file to stdout.

Designed for use with a BEDTools pipeline.
"""

import optparse
import bedparser
import sys
import os

op = optparse.OptionParser(usage='asdf')
op.add_option('-i', type=str, dest='input', help='input bed file, if "stdin" will accept input from stdin.')
op.add_option('-f', type=int, dest='flank', help='flanking region to take on each side')
op.add_option('-l', type=int, dest='left', help='specify left side flanking region')
op.add_option('-r', type=int, dest='right', help='specify right side flanking region')
op.add_option('--buffer', type=int, default=0, dest='buffer', 
               help='part of flanking region closest to feature that you want to exclude from flanking region')
__doc__ += op.format_help()

def main(options):
    if (options.left or options.right) and options.flank:
        raise ValueError, '-f cannot be specified if -r or -l is used.'
        sys.exit(1)

    if options.flank:
        options.left = options.flank
        options.right = options.flank

    if not options.input:
        raise ValueError, 'Need input file.'

    if options.input == 'stdin':
        options.input = sys.stdin

    if not options.left or not options.right:
        raise ValueError, 'Flanking region not specified'
        

    s = '%s\t%s\t%s'
    for i in bedparser.bedfile(options.input):
        print s % (i.chr, i.start - options.left - options.buffer, i.start - options.buffer)
        print s % (i.chr, i.stop + options.buffer, i.stop + options.right + options.buffer)

if __name__ == "__main__":
    
    options,args = op.parse_args()
    main(options)
