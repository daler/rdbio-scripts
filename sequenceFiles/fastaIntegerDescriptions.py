#!/usr/bin/python
__doc__="""
Script to replace sequence descriptions in a FASTA file with an integer. Useful
for when the Genome Browser packs a lot of info into a fasta record description
but that info crashes downstream programs.

"""

import optparse
import sys 
from cStringIO import StringIO
from Bio import SeqIO

op = optparse.OptionParser(usage='')
op.add_option('-i', dest='infn', help='Input FASTA file.')
op.add_option('-o', dest='outfn', help='Output FASTA file, with integers as descriptions')
__doc__ += op.format_help()

def fastaIntegerDescriptions(infile,outfile):
    def seq_iterator(handle):
        count = 0
        for i in SeqIO.parse(open(fin),'fasta'):
            c = str(count)
            i.description = c
            i.id = c
            i.name = c
            count += 1
            yield i
    SeqIO.write(seq_iterator(handle), outfile, 'fasta')


if __name__ == "__main__":
    options,args = op.parse_args()
    if options.infn is None:
        infile = sys.stdin
    else:
        infile = open(options.infn)
    if options.outfn is None:
        outfile = sys.stdout
    else:
        outfile = open(options.outfn, 'w')
    fastaIntegerDescriptions(infile,outfile)
    
    if options.outfn is not None:
        outfile.close()
