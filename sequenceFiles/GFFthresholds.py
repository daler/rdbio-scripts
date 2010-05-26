#!/usr/bin/python
"""
This script applies the following criteria to a GFF file: 4 probes separated by
a max of 500 bp that are greater than 2.5sd greater than the mean log ratio.

(parameters are what Ren et al used)

Input is a Nimblegen GFF file.
Output is a BED file.
Test data are found in tests/GFFthresholds
"""
from pylab import csv2rec
import logging
import sys

def nodate(x):
    return str(x)

def GFFthreshold(infn,outbed):
    """
    Thresholds the values in the GFF file *infn* and exports
    the results to the BED file *outbed*.
    """
    converterd = {'probe':nodate,'a':nodate,'b':nodate}
    logging.debug('reading GFF into record array')
    a = csv2rec(infn, 
                delimiter='\t', 
                names=('chr','prog','id','start','stop','ratio','a','b','probe'),
                converterd=converterd)
    logging.debug('sorting record array')
    a.sort(order=('chr','start'))
    fout = open(outbed,'w')
    m = a.ratio.mean()
    std = a.ratio.std()
    thresh = m + 2.5 * std
    allregions = []
    region = []
    lastchr = a.chr[0]
    lastpos = None
    count = 0

    for data in a:
        if data.ratio < thresh:
            continue

        if lastpos is None:
            dist = 0
        else:
            dist = data.start - lastpos
        
        logging.debug('region is currently')
        for i in region:
            logging.debug('\t%s' % i)
        logging.debug('this data: %s' % data)
        logging.debug('dist from last: %s' % dist)
    
        if dist > 500 or data.chr != lastchr:
            
            logging.debug('\ndist > 500; checking region len')
            logging.debug('regionlen: %s' % len(region))
            for i in region:
                logging.debug('\t%s' % i )
            if len(region) < 4:
                logging.debug('region not long enough, erasing')
            else:
                logging.debug('region is long enough!!!!')
                logging.debug('region to be exported is')
                for i in region:
                    logging.debug('\t%s' % i)
                chr = region[0].chr
                start = region[0].start
                stop = region[-1].stop
                fout.write('%s\t%s\t%s\n' % (chr,start,stop))
                count += 1
            region = []

        lastpos = data.stop
        lastchr = data.chr
        logging.debug('adding %s to region' % data)
        region.append(data)

    if len(region) >= 4:
        logging.debug('last region will be exported')
        logging.debug('region to be exported is')
        for i in region:
            logging.debug('\t%s' % i)
        
        chr = region[0].chr
        start = region[0].start
        stop = region[-1].stop
        fout.write('%s\t%s\t%s\n' % (chr,start,stop))
        count += 1

    else:
        logging.debug('last region not long enough')

    fout.close()
    logging.debug('Number of enriched regions: %s' % count)
    logging.debug('using threshold: %s' % thresh)


if __name__ == "__main__":
    
    import optparse
    op = optparse.OptionParser()
    op.add_option('-i',dest='infn', help='GFF file input')
    op.add_option('-o',dest='outbed',help='Output bed file')
    op.add_option('--test',dest='test',action='store_true',help='run tests')
    op.add_option('--verbose',dest='verbose',action='store_true', help='print debugging info')
    options,args = op.parse_args()

    reqargs = ['infn','outbed']
    fail = False
    for rq in reqargs:
        if not getattr(options,rq):
            if not fail:
                print 
            print '\tMissing required arg %s' % rq.upper()
            fail = True
    if fail:
        print '\n'
        op.print_help()
        sys.exit()

    if options.test:
        options.infn = 'tests/GFFthresholds/test.gff'
        options.outbed = 'tests/GFFthresholds/test.gff.out.bed'
    if options.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)

    thresh = GFFthreshold(options.infn, options.outbed)


