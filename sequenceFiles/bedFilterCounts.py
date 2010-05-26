#!/usr/bin/python
__doc__ = """Takes as input the output of ``intersectBed`` with the ``-c``
flag.  Returns a new bed file with only those records that have a count (the
last field) > N.

"""
import optparse, sys
usage = ''
op = optparse.OptionParser(usage='')
op.add_option('-i', dest='infn', help='Input bed file, can be "stdin"')
op.add_option('-n', dest='n', help='Return only features with counts greater than N.')
op.add_option('-o', dest='outfn', help='If not specified, defaults to stdout')
op.usage = usage
__doc__ += op.format_help()


def bedFilterCounts(infn, n, outfn=None):
    """
    *infn*
        
        Already opened file-like object

    *n*

        Integer, only features that have a count *greater* than *n* will be 
        reported.

    *outfn*

        Already opened file-like object for writing.
        """
    try:
        n = int(n)
    except ValueError:
        print '-n must be an integer'
        sys.exit()
    if outfn is None:
        fout = sys.stdout
    else:
        fout = open(outfn, 'w')
    if infn == 'stdin':
        fin = sys.stdin
    else:
        fin = open(infn)
    
    for line in fin:
        L = line.rstrip().split('\t')
        try:
            count = int(L[-1])
        except ValueError:
            print "The last column doesn't appear to be an integer"
            sys.exit(1)
        if count > n:
            fout.write(line)

    if outfn is not None:
        fout.close()

if __name__ == "__main__":
    options, args = op.parse_args()
    bedFilterCounts(options.infn, options.n, options.outfn)
