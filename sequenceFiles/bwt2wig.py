"""
Script to parse a Bowtie output file into a wig file, where the value at each
position in the wig file corresponds to the number of [overlapping] reads
mapping to that location. 

In this implementation, writing the WIG files is pretty slow -- on the order of
3 mins for a single chromosome.

Some ways to optimize:

    * There are a lot of zeros.  Should probably find the non-zero indexes and
      only write those out to file. Custom iterator prob the way to go for this.
    * Only write to file when you really mean it . . . so batch together a bunch 
      of [non-zero] values into a string, then write that long string to file.

Created: Jan 2010 Ryan Dale
"""
import optparse
from numpy import zeros
import sys
import time
usage = """

Script to parse bowtie-format files into WIG-format files.  

Values of each base in the wig file correspond to the number of stacked reads
at that position.  If --chrom is not specified, all chromosomes will be used.

Separate files will be created for each chromosome."""
op = optparse.OptionParser(usage=usage)

op.add_option('-i', dest='input', help='Input Bowtie-format file')
op.add_option('--prefix', dest='prefix', help='Prefix for output WIG filenames.  '
                                              'These filenames will also contain '
                                              'the chromosome name (e.g., PREFIX.chr3L.wig)')
op.add_option('--chrom', dest='chrom', help='Optional chromosome name to use.  '
                                            'If specified, will only output a '
                                            'single file for this chromosome and '
                                            'will ignore the others.')
options,args = op.parse_args()
fin = open(options.input)

# got these values from the chromInfo table on UCSC
chromLimits = {
'chr2L'      : 23011544,
'chr2LHet'   : 368872,
'chr2R'      : 21146708,
'chr2RHet'   : 3288761,
'chr3L'      : 24543557,
'chr3LHet'   : 2555491,
'chr3R'      : 27905053,
'chr3RHet'   : 2517507,
'chr4'       : 1351857,
'chrU'       : 10049037,
'chrUextra'  : 29004656,
'chrX'       : 22422827,
'chrXHet'    : 204112,
'chrYHet'    : 347038,
'chrM'       : 19517
}

# prepare dictionary of initialized arrays for each chrom
chromDict = {}
for chrom, chrom_len in chromLimits.items():
    chromDict[chrom] = zeros( (chrom_len,), dtype=int)  # int dtype to save on memory

# Parse the bowtie output file
t0 = time.time()
counter = 0
for line in fin:
    counter += 1
    if counter % 10000 == 0:
        print '\r%s reads processed'%counter,
        sys.stdout.flush()
    L = line.split('\t')
    strand = L[1]
    chrom = L[2]
    start = int(L[3])
    seq = L[4]
    stop = start + len(seq)-1

    # simply increment the values at positions that overlap this read by 1
    chromDict[chrom][start:stop] += 1 

fin.close()
t1 = time.time()
print '\nParsing Bowtie file took %ds.' % (t1-t0)


def iterator(x):
    """
    Iterator that yields bunches of non-zero lines to be written to a wig file.
    *x* is a NumPy array.
    """
    THRESH = 100
    s = []
    start = 0
    run_of_zeros = 0
    last_i_was_zero = False
    for ind,i in enumerate(x):
        if i == 0:
            if last_i_was_zero:
                run_of_zeros += 1
            last_i_was_zero = True
        else:
            last_i_was_zero = False
        s.append(str(i))
        if run_of_zeros == THRESH:
            if run_of_zeros < len(s):
                subset = s[:-THRESH+1]
                if subset > 1:
                    yield 'stepSize=1 start=%s\n' % (start) + '\n'.join(subset)+'\n'
            s = []
            last_i_was_zero = False
            run_of_zeros = 0
            start = ind
            


# Write each chromosome's data to a separate file (a la MACS)
print '\nWriting wig file...'
for chrom,values in chromDict.iteritems():
    t0 = time.time()
    if options.chrom is not None:
        if chrom != options.chrom:
            continue
    print chrom,
    fout = open('%s.%s.wig'%(options.prefix,chrom),'w')
    sys.stdout.flush()
    fout.write('track type=wiggle_0 name=%s_coverage\n' % chrom)
    fout.write('fixedStep chrom=%s start=1 step=1\n' % chrom)
    for chunk in iterator(values):
        fout.write(chunk)
    t1 = time.time()
    fout.close()
    print '(%ds)'%(t1-t0)
