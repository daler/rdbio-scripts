#!/usr/bin/python

__doc__ = """
Script to split a FASTA file with (possibly) more than one line per record in
to a file with all bases on one line.  New files are created for each record in
the fasta file.  Supply a --prefix=DEST argument to save the files in a
directory other than the current directory.

This format is needed for, e.g., BioProspector.

Input can be file or stdin; output is stdout.
"""

import sys
import optparse
import os
from Bio import SeqIO

op = optparse.OptionParser()
op.add_option('-i', dest='infn', help='Input FASTA file.')
op.add_option('--prefix', dest='prefix', default='', help='Directory in which to store split files.')
op.add_option('--verbose', dest='verbose',action='store_true',help='Verbose mode')
__doc__ += op.format_help()

def main(options):
    handle = open(options.infn)
    extension = os.path.splitext(options.infn)[-1] 
    
    for record in SeqIO.parse(handle, 'fasta'):
        fout = os.path.join(options.prefix,record.description+extension)
        fout = open(fout,'w')
        seq = record.seq.tostring()
        description = record.description
        if options.verbose:
            print description
        fout.write('>%s\n' % description)
        fout.write(seq)
        fout.close()

if __name__ == "__main__":
    options,args = op.parse_args()
    main(options)
