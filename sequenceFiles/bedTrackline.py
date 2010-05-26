#!/usr/bin/python

__doc__ = """Prepends a track line to the bed file. Can be piped to (use -i
stdin) or from (don't specify -o)."""

import optparse,sys

op = optparse.OptionParser()
op.add_option('-i', dest='infn', help='Input bed file, can be "stdin"')
op.add_option('--trackline', dest='trackline', help='Replaces (if exists) or adds new track line')
op.add_option('-o', dest='outfn', help='If not specified, defaults to stdout')
__doc__ += op.format_help()

def bedTrackline(infn, trackline, outfn=None):
    """Adds the track name *trackline* to the beginning of 
    *infn*, outputting to *outfn* (which, if None, defaults to stdout).
    
    *infn* can be "stdin".
    """ 
    if outfn is None:
        fout = sys.stdout
    else:
        fout = open(outfn, 'w')
    if infn == 'stdin':
        fin = sys.stdin
    else:
        fin = open(infn)
    fout.write('%s\n' % trackline)
    s = fin.readline()
    if not s.startswith('track'):
        fout.write(s)
    s = fin.read()
    fout.write(s)

    if outfn is not None:
        fout.close()
    if infn is not 'stdin':
        fin.close()

if __name__ == "__main__":
    options, args = op.parse_args()
    bedTrackline(options.infn, options.trackline, options.outfn)
