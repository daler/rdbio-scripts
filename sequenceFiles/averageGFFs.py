#!/usr/bin/python

# Ryan Dale
# Created June 2009

import optparse, sys, logging
from numpy import median, abs, sum, array, mean

description = """This script averages the data among several GFF files of the
same cell type and results in a single GFF file."""
usage = "\n\n%prog file1 [file2 fileN] [options] -o OUTPUTFILE"
op = optparse.OptionParser(usage=usage,description=description)
op.add_option('-o', 
              help='output file', 
              dest='outputfile')
op.add_option('--biweight-scaling', 
    help='''subtract biweight M-estimator 
from the averaged values (default is %default)''',
    dest='biweight',
    action='store_true',
    default=True)
op.add_option('-c', 
    help='parameter "c" for biweight scaling (default is %default)', 
    dest='c',
    default=5, 
    type='float')
op.add_option('--epsilon',
    help='epsilon value for biweight scaling (default is %default)',
    dest='epsilon',
    type='float',
    default=1e-4)
options,args = op.parse_args()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('log')

if len(args) < 1: 
    logger.warning('''need at least one input (and preferably multiple) 
files to work on! (use -h switch to view help)''')
    sys.exit()

try:
    options.outputfile
except AttributeError:
    logger.warning('No output file specified! (use -h switch to view help)')
    sys.exit()


# These two dictionaries will be keyed by probe ID
logratios = {}
otherstuff = {}

def biweight(x,c=5.0,epsilon=1e-5):
    '''Computes the Tukey biweight M-estimator.'''
    x = array(x)
    m = median(x)
    s = median(abs(x-m))
    u = (x-m) / (c*s + epsilon)
    w = x*0
    i = abs(u) <= 1
    w[i] = ((1-u**2)**2)[i]
    t = sum(w*x) / sum(w)
    return t

for inputfile in args:
    logger.info('Reading in %s' % inputfile)
    f = open(inputfile)
    for line in f:
        if line.startswith('#'):
            continue
        L = line.strip().split('\t')
        key = L[-1]
        log2ratio = float(L[5])
        vals = logratios.setdefault(key,[])
        vals.append(log2ratio)
        try:
            otherstuff[key]
        except KeyError:
            otherstuff[key] = L

# get the mean for each item
logger.info('Getting mean for each feature')
for key,value in logratios.iteritems():
    logratios[key] = mean(value)

# calculate the biweight scaling if requested

if options.biweight:
    logging.info('Calculating biweight scaling')
    global_biweight = biweight(logratios.values())

else:
    logger.info('No biweight scaling requested, skipping')

# construct lines for output file.
logger.info('creating %s' % options.outputfile)
fout = open(options.outputfile, 'w')
for key,value in logratios.iteritems():
    linelist = otherstuff[key]
    if options.biweight:
        linelist[5] = str(value-global_biweight)
    else:
        linelist[5] = str(value)
    line = '\t'.join(linelist) + '\n'
    fout.write(line)

fout.close()
logger.info('DONE!')







