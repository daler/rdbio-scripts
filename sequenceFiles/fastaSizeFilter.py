#!/usr/bin/python
__doc__ = """
Filters a FASTA file by sequence length.  Only sequences > L will be sent to the 
output file.

``--hist`` will plot a histogram of the lengths.

Prints to stdout.

"""

import optparse
import sys
from Bio import SeqIO
from cStringIO import StringIO
from matplotlib import pyplot as p

op = optparse.OptionParser()
op.add_option('-i', dest='infn', help='Input FASTA-format file. If unspecified, use stdin.')
op.add_option('-o', dest='outfn', help='Output FASTA-format file, with records of length *L* or greater. If unspecified, use stdout.')
op.add_option('--hist', action='store_true', dest='hist',help='Plot histogram of pre-filtered records')
op.add_option('--hist-kwargs', dest='histkwargs', help='Comma-separated keyword args to pass to matplotlib histogram function')
op.add_option('-L', type=int, dest='L', help='Seqs < L will be filtered out')
op.add_option('--upper', action='store_true', dest='upper',help='Force letters to be uppercase')

if not options.input:
    raise ValueError, 'Need an input file'
if not options.L:
    raise ValueError, 'Need a length to filter by'

def _length_filter(handle,L):
    """Iterator that yields sequences >= L."""
    for i in SeqIO.parse(handle,'fasta'):
        if len(i) >= L:
            s = i.seq.tostring()
            s = s.upper()
            seq = i.seq.tomutable()
            seq[:] = s
            i.seq = seq       
            yield i

def _fasta_lengths(handle):
    """Returns a list of seq lengths"""
    lengths = []
    for i in SeqIO.parse(open(fin), 'fasta'):
        lengths.append(len(i))
    return lengths

def _parse_kwargs(kw):
    '''Return a dictionary parsed from a string of comma-separated
    kwargs'''
    d = {}
    try:
        L = kw.split(',')
        for i in L:
            key,value = i.split('=')
            try: 
                value = float(value)
            except ValueError:
                value = value
            d[key] = value

    except AttributeError:  # NoneType has no split()
        pass
    return d

def fastaSizeFilter(infile, outfile, L, hist=False, histkwargs=None, upper=False):
    """
    Filter *infile* so that only records of length *L* or greate are reported
    in *outfile*.

    *infile*

        Already-open input FASTA file

    *outfile*

        Already-open output file.

    *L*

        Length, in bases, of the minimum length a record should be.  Only
        records greater than this will be reported

    *hist*

        Default False.  If True, will plot a histogram of the lengths of records in 
        in *infile*.
    """
    seq_iterator = _length_filter(infile, L)
    SeqIO.write(seq_iterator, outfile, 'fasta')
    kwargs = _parse_kwargs(histkwargs)
    if hist:
        fig = p.figure()
        ax = fig.add_subplot(111)
        ax.hist(_fasta_lengths(handle),**kwargs)
        ax.set_ylabel('count')
        ax.set_xlabel('Sequence length')
        p.show()

if __name__ == "__main__":
    options,args = op.parse_args()
    if options.infn is None:
        infile = sys.stdin
    else:
        infile = open(options.infn)
    fastaSizeFilter(infile, outfile, options.L, options.hist, options.histkwargs, options.upper)
