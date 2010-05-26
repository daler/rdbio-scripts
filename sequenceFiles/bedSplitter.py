#!/usr/bin/python
import optparse, sys, bedparser, os, tempfile
__doc__ = """
This script splits a BED-format file into separate tracks
"""

op = optparse.OptionParser(usage='')
op.add_option('-i', dest='infn', help='Input file.  stdin support when bedparser gets it.')
op.add_option('-o', dest='outfn', help='Output file.  If unspecified, use stdout.')
__doc__ += op.format_help()


def bedSplitter(infile, outfile):
    """
    Bins BED features with common names into individual tracks.

    *infile*

        An open file-like object

    *outfile*

        A file-like object open for writing.
    """
    # first, make sure it's sorted by name.
    outfile = open(outfile,'w')
    tmp = tempfile.mktemp()
    cmd = 'sort -k 4 %s > %s' % (infile, tmp)
    print cmd
    os.system(cmd)
    trackname = None
    lastname = None
    for i in bedparser.bedfile(tmp):
        if i.name != lastname:
            trackname = i.name
            outfile.write('track name=%s description=%s itemRgb=1\n' % (i.name,i.name))
        outfile.write(i.tostring())
        lastname = i.name


if __name__ == "__main__":
    options,args = op.parse_args()

    bedSplitter(options.infn, options.outfn)


