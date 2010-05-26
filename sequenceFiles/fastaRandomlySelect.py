#!/usr/bin/python

__doc__ = """Randomly selects *n* fasta sequences.

"""
import sys
import random
import optparse
from Bio import SeqIO
from cStringIO import StringIO

op = optparse.OptionParser(usage='')
op.add_option('-i',dest='infn',help='Input FASTA file.  If unspecified, reads from stdin.')
op.add_option('-o',dest='outfn',help='Output FASTA file.  If unspecified, writes to stdout.')
op.add_option('-n',dest='n', type=int, help='Number of sequences to return.')
__doc__ += op.format_help()

def fastaRandomlySelect(infile,outfile,n):
    """
    *infn*

        Input FASTA file (an open file-like object).  If None, reads
        from stdin.

    *n*

        Integer, number of sequences to select.

    *outfn*

        Output FASTA file containing *n* sequences.  (an open file-like object)
        If None, prints to stdout.
    """
    # strategy: read in all as a list, randomly shuffle, 
    # take the first N.

    seqs = SeqIO.parse(infile,'fasta')
    seqs = list(seqs)
    random.shuffle(seqs)
    SeqIO.write(seqs[:options.n], outfile, 'fasta')

if __name__ == "__main__":
    options,args = op.parse_args()
    if options.infn is None:
        infile = sys.stdin
    else:
        infile = open(options.infn)
    if options.outfn is None:
        outfile = sys.stdout
    else:
        outfile = open(options.outfn,'w')
    if not options.n:
        print 'Need to specify -n'
        sys.exit()

    fastaRandomlySelect(infile,outfile,options.n)

    if options.outfn is not None:
        outfile.close()
