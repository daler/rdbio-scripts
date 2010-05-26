#!/usr/bin/python

__doc__ = """
Filters a BED file by value, given the ``--min`` and ``--max``
argmuments.  Does not yet accept stdin (need to get bedparser
to do that).

Writes to stdout."""

import sys
import optparse
import bedparser

op = optparse.OptionParser(usage='%prog [options]')
op.add_option('-i',dest='infn',help='Input BED file.')
op.add_option('--min',type=float,dest='min',help='Minimum value to accept.')
op.add_option('--max',type=float,dest='max',help='Maximum value to accept.')
op.add_option('--bedgraph', action='store_true', dest='bedgraph', help='Work on a bedGraph file instead of bed')
op.add_option('-o',dest='outfn',help='output filtered BED file')
__doc__ += op.format_help()

def bedValueFilter(fn, minval, maxval, outfn, bedgraph=False):
    """Filters a bed file by value.  *fn* is a filename, while *outfn* is a
    file-like handle.""" 
    if maxval is None:
        maxval = 1e30
    if minval is None:
        minval = -1e30
    
    if bedgraph:
        iterator = bedparser.bedgraph(fn)
    else:
        iterator = bedparser.bedfile(fn)

    for i in iterator:
        if i.value is None:
            print "No value for this feature (%s)" % i
            sys.exit(1)
        if minval < i.value < maxval:
            fout.write('%s\t%s\t%s\t%s\n' % (i.chr, i.start,i.stop,i.value))
    if outfn is not None:
        fout.close()

if __name__ == "__main__":
    options,args = op.parse_args()
    if options.outfn is None:
        fout = sys.stdout
    else:
        fout = open(options.outfn,'w')

    bedValueFilter(options.infn, options.min, options.max, fout,bedgraph=options.bedgraph)

