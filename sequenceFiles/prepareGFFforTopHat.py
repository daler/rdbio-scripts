"""
test_chr test gene 101 300 . + . ID=test1gene
test_chr test mRNA 101 300 . + . ID=test1mRNA;Parent=test1gene
test_chr test exon 101 148 . + . ID=test1exon1;Parent=test1mRNA
test_chr test exon 253 300 . + . ID=test1exon2;Parent=test1mRNA
"""
import optparse
from textwrap import dedent as dd
import logging
import subprocess

logging.basicConfig(level=logging.DEBUG)

op = optparse.OptionParser()
op.add_option('-i', dest='infile', help='Input GFF file.')
op.add_option('-o', dest='outfile',
              help='Output GFF file, suitable for use with TopHat')
op.add_option('--index', 
              dest='index',
              help=dd('''
                   index name for Bowtie index. 
                   bowtie-inspect and the index 
                   must be on your path.'''))
op.add_option('--no-auto-chr',dest='noautochr',
              default=False,
              action='store_true',
              help="""Disable automatic prepending of "chr" to output file
              chromosome names if not already present.""")
options,args = op.parse_args()

p = subprocess.Popen('bowtie-inspect --names %s' % options.index, 
                     shell=True, stdout=subprocess.PIPE)
chroms = p.communicate()[0]
chroms = chroms.split('\n')
chroms = [i for i in chroms if len(i) > 0]
logging.debug('chromosomes found in bowtie index %s were \n\t%s' \
                % (options.index, '\n\t'.join(chroms)))


def parsedescription(description):
    '''parses the description field of a GFF3 file; returns a string of
    only the ID and Parent info'''
    d = {}
    items = description.split(';')
    for i in items:
        field,data = i.split('=')
        d[field]=data
    try:
        s ='ID=%s;Parent=%s' % (d['ID'],d['Parent'])
    except KeyError:
        s = 'ID=%s' % d['ID']
    return s

featuretypes =  ['gene','mRNA','exon','ncRNA','miRNA','rRNA','snoRNA','snRNA','tRNA']
count = 0
outfile = open(options.outfile,'w')
for line in open(options.infile):
    if line.startswith('#'):
        continue # comments
    if line.startswith(">"):
        break # you've reached the FASTA part of the file.
    L = line.rstrip().split('\t')

    chrom,name,featuretype,start,stop,value,strand,_,description = L
    
    if 'chr' not in chrom and not options.noautochr:
        chrom = 'chr%s' % chrom
    
    if chrom not in chroms: # wasn't in the bowtie index, so skip.
        continue
    if featuretype not in featuretypes:
        continue
    parsed_desc = parsedescription(description)

    if 'RNA' in featuretype:
        featuretype = 'mRNA'

    outline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom,
                                                        name,
                                                        featuretype,
                                                        start,
                                                        stop,
                                                        value,
                                                        strand,
                                                        _,
                                                        parsed_desc)
    outfile.write(outline)
    count += 1
logging.info('Wrote %s features to %s' % (count, options.outfile))

     

    
