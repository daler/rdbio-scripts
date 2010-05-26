#!/usr/bin/python

__doc__ = """
Filters a BED file by size, given the ``--min`` and ``--max``
argmuments."""

import sys
import optparse
import bedparser

op = optparse.OptionParser(usage='%prog [options]')
op.add_option('-i',dest='infn',help='Input BED file')
op.add_option('-o', dest='outfn', help='Output BED file; if unspecified will use stdout')
op.add_option('--min',type=float,dest='min',help='Minimum length to accept.')
op.add_option('--max',type=float,dest='max',help='Maximum length to accept.')
__doc__ += op.format_help()

def bedSizeFilter(infile, outfile, minlen=None, maxlen=None):
    if maxlen is None:
        maxlen = 1e30
    if minlen is None:
        minlen = -1e30
    for i in bedparser.bedfile(infile):
        length = abs(i.start-i.stop)
        if minlen < length < maxlen:
            outfile.write('%s\t%s\t%s\n' % (i.chr, i.start,i.stop))

if __name__ == "__main__":
    options,args = op.parse_args()
    if options.outfn is None:
        outfile = sys.stdout
    else:
        outfile = open(options.outfn,'w')
    if options.infn is None:
        infile = sys.stdin
    else:
        infile = options.infn
    bedSizeFilter(infile, outfile, options.min, options.max)
    if options.outfn is not None:
        outfile.close()
