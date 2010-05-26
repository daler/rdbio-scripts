#!/usr/bin/python

__doc__ = """Script to replace any lowercase values with an N.  Prints to
stdout.  Useful for when you need RepeatMasked sequences from a fasta file
downloaded from UCSC.

"""

import sys
import optparse
from cStringIO import StringIO
from Bio import SeqIO
from string import lowercase

op = optparse.OptionParser(usage='')
op.add_option('-i', dest='infn', help='Input FASTA file.  Can be "stdin".')
op.add_option('-o',dest='outfn',help='Output FASTA file, with lowercase replaced with "N"')

__doc__ += op.format_help()


def fastaLowercaseToNs(infn,outfn):

    if infn == 'stdin':
        handle = sys.stdin
    else:
        handle = open(infn)

    if outfn is None:
        fout = sys.stdout
    else:
        fout = open(outfn,'w')


    def N_seq_iterator(handle):
        """Generator function to return an iterator.  Pass this 
        to SeqIO.write, which takes an iterator as an argument.
        
        Returns a SeqRecord object with lowercase letters replaced
        with N."""
        for i in SeqIO.parse(handle, 'fasta'):
            s = i.seq.tostring()
            new_s = []
            for b in s:
                if b in lowercase:
                    new_s.append('N')
                else:
                    new_s.append(b)

            new_s = ''.join(new_s)

            seq = i.seq.tomutable()
            seq[:] = new_s
            i.seq = seq
            yield i

    seqs = N_seq_iterator(handle)
    SeqIO.write(seqs, fout, 'fasta')
    sys.stdout.flush()

    if outfn is not None:
        fout.close()

if __name__ == "__main__":
    options,args = op.parse_args()
    fastaLowercaseToNs(options.infn,options.outfn)
