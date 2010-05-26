#!/usr/bin/python

__doc__ = """
Script to convert FASTA files with (possibly) more than one line per record in
to a file with all bases/amino acids on one line.

This format is needed for, e.g., BioProspector.

Input can be file or stdin; output is stdout.

Requires BioPython.
"""

import sys
import optparse
from Bio import SeqIO

op = optparse.OptionParser()
op.add_option('-i', dest='infn', help='Input FASTA file.  Can be "stdin".')
__doc__ += op.format_help()

def main(options):
    if options.infn == 'stdin':
        handle = sys.stdin
    else:
        handle = open(options.infn)
    fout = sys.stdout
    for record in SeqIO.parse(handle, 'fasta'):
        seq = record.seq.tostring()
        description = record.description
        fout.write('>%s\n' % description)
        fout.write(seq+'\n')

if __name__ == "__main__":
    options,args = op.parse_args()
    main(options)
