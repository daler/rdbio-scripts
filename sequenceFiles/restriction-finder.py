#!/usr/bin/python
usage ="""

Makes a BED file of the cut sites of a specified restriction enzyme.

Example usage:
    
    # Get BED file of DpnI sites in dm3.fa
    python restriction-finder.py --fasta dm3.fa --enzyme DpnI --bed DpnI-sites.bed

    # can pipe to BedTools to get, e.g, sites in genes::
    python restriction-finder.py --fasta myfasta.fa --enzyme DpnI | intersectBed -a stdin -b genes.bed > DpnI-in-genes.bed


Created 13 Aug 2010 by Ryan Dale"""
try:
    from Bio import SeqIO
    from Bio import Restriction
except ImportError:
    sys.stderr.write("\nPlease install BioPython to use this script <http://biopython.org/wiki/Biopython>\n")
import optparse
import sys
import os

op = optparse.OptionParser(usage=usage)
op.add_option('--fasta', help='Required FASTA file containing sequences to search')
op.add_option('--enzyme', help='Required enzyme name, case sensitive (e.g., DpnI or EcoRI)')
op.add_option('--bed',help='BED file to create. If not specified, output will print to stdout.')
options,args = op.parse_args()

# Input error checking...
def err(s):
    op.print_help()
    sys.stderr.write('\n***ERROR: %s***\n\n'%s)
    sys.exit(1)
# Hack to import just the enzyme you want from the Restriction module
if options.enzyme is None:
    err('Please specify an enzyme with --enzyme')
if options.fasta is None:
    err('Please specify a FASTA file with --fasta')
try:
    exec('from Bio.Restriction import %s as restr' % options.enzyme)
except ImportError:
    err('No restriction enzyme "%s" found!' % options.enzyme)

if not os.path.exists(options.fasta):
    err('FASTA file %s not found'%options.fasta)

if options.bed is None:
    fout = sys.stdout
else:
    fout = open(options.bed,'w')


# Let BioPython do the work...
parser = SeqIO.parse(options.fasta,'fasta')
for chrom in parser:
    sys.stderr.write(chrom.name+'\n')
    hits = restr.search(chrom.seq)
    for hit in hits:
        values = [chrom.name,
                  str(hit),
                  str(hit+1)]
        fout.write('\t'.join(values)+'\n')
        fout.flush()
if options.bed is not None:
    fout.close()
