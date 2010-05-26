#!/usr/bin/python
"""
Script to convert SAM, BED or Bowtie-formatted files into a "pileup" or
"coverage" wig file.

***IMPORTANT***

    Assumes that the input file is sorted by chromosome, then by start
    position.  For a Bowtie-format file, sort your file with::

        sort -k 3,3 -k 4n INPUTFILE > OUTPUTFILE

    For a BED file, use::

        sort -k 1,1 -k 2n INPUTFILE > OUTPUTFILE

How does it work?
-----------------
The strategy is to bin together contiguously overlapping features into clusters
of reads, computing the cumulative coverage on each cluster.  Then the cluster
is written to file as WIG format and removed from memory. The result is that
only one cluster is in memory at a time.  

In contrast to previous implementations that intitialized an array the size of
the genome, long stretches of empty genome are simply not checked here -- only
the features in the input file are used.

Singleton reads not showing up?
-------------------------------

Since only non-zero values are written to the WIG file, the auto-scaling of the
UCSC browser means that the lowest y-value visible in the track will be 1.  As
a result, the singleton reads will appear to not have any wiggle track (since
the y-axis starts at 1).  To force them to show up, you can manually add a zero
value to the beginning of the WIG file, after the track line, like so::

    fixedStep chrom=chr4 start=1
    0

Then you can manually set the y-axis scaling in the browser to include zero,
which will make the singleton reads show up.

Other formats
-------------
Support for other input file types can be implemented by writing another
iterator that returns chrom, start, stop (start and stop must be integers) and
adding that to the dispatch dictionary.  


Created Jan 2010 Ryan Dale
Updated Mar 8 2010 Ryan Dale added support for SAM format
Updated Apr 15 2010 Ryan Dale added support for scale factor
"""

usage = """
Convert a pre-sorted SAM, BED or Bowtie file into WIG format.

***IMPORTANT***

    Assumes that the input file is sorted by chromosome, then by start
    position.  For a Bowtie-format file, sort your file with:

        sort -k 3,3 -k 4n INPUTFILE > OUTPUTFILE

    For a BED file, use::

        sort -k 1,1 -k 2n INPUTFILE > OUTPUTFILE
"""

import optparse
import time
import sys

op = optparse.OptionParser(usage=usage)
op.add_option('-i',     dest='input',  help='Input file to be converted to WIG format')
op.add_option('--type', dest='type',   help='Type of input file, one of "bed", "sam" or "bowtie". Required.')
op.add_option('-o',     dest='output', help='Optional output file to write to. If unspecified, writes to stdout.')
op.add_option('--scale',dest='scale',  help='Optional scale factor to multiply all WIG values by. Useful '
                                            'for comparing datasets with differing library sizes', type=float,
                                            default=1.0)

# NOTE: Support for other input file types can be implemented by writing
# another iterator that returns chrom, start, stop (start and stop must be
# integers) and adding that to the dispatch dictionary.

def samfile_iterator(fn):
    """
    Iterates through a SAM file *fn*, returning chrom,start,stop for each
    line where there was a mapped read.  start and stop are converted to ints.

    Ignores unmapped reads.
    """
    for line in open(fn):
        L = line.strip().split('\t')
        flag = int(L[1])
        if flag & 0x0004:
            continue
        
        chrom = L[2]
        start = int(L[3])-1
        stop = start + len(L[9])
        yield (chrom,start,stop)
            

def bedfile_iterator(fn):
    """
    Iterates through a BED file *fn*, returning chrom,start,stop for each line
    (start and stop are converted to int)

    Track lines are gracefully ignored.
    """
    for line in open(fn):
        L = line.strip().split()
        chrom = L[0]
        try:
            start = int(L[1])
        except ValueError:
            if line.startswith('track'):
                continue
            if line.startswith('browser'):
                continue
            else:
                raise ValueError,'mis-formatted bed line'
        stop = int(L[2])
        yield (chrom,start,stop)

def bowtie_iterator(fn):
    """
    Iterates through a Bowtie-format file *fn*, returning chrom,start,stop for
    each line (start and stop are converted to int)
    """
    for line in open(fn):
        L = line.strip().split()
        chrom = L[2]
        start = int(L[3])
        stop = start + len(L[4])-1
        yield (chrom,start,stop)

# Makes the connection between filetype name and the iterator to use for that
# filetype
dispatch_dict = {
                 'bed':bedfile_iterator,
                 'bowtie':bowtie_iterator,
                 'sam':samfile_iterator,
                }

def clusters(infn, filetype):
    """
    Yields clusters of overlapping reads along with the chromsome and cluster boundaries.

    Return value is of the form (chrom, cluster_start, cluster_stop, features)
    where *features* is a list of (start,stop) tuples.
    """

    # Decide which iterator to use for this filetype.
    try:
        iterator = dispatch_dict[filetype](infn)
    except KeyError:
        raise ValueError, "Input filetype %s not supported; only %s currently supported" % (filetype,dispatch_dict.keys())

    last_chrom = None
    for chrom,start,stop in iterator:

        # If we're on a new chromosome, reset everything.
        if chrom != last_chrom:
            sys.stderr.write('%s\n'%chrom)
            features = []
            cluster_start = start
            cluster_stop = stop
        last_chrom = chrom 

        # Check for overlap with the current cluster.  If it overlaps, then
        # extend the cluster limits as needed, and add the features to the list
        # of overlapping features (this list is the cluster itself)
        if start <= cluster_stop:
            if start <= cluster_start:
                cluster_start = start
            if stop >= cluster_stop:
                cluster_stop = stop
            features.append( (start,stop) )

        else:
            # This feature doesn't overlap, so return the existing cluster...
            yield chrom,cluster_start,cluster_stop,features

            # ...and start a new cluster with this feature.
            cluster_start = start
            cluster_stop = stop
            features = [(start,stop)]
            

    # Return the last cluster (since it won't be yielded by the else-clause above)
    yield chrom,cluster_start, cluster_stop, features

def make_wig(infn, filetype, outfn=None):
    """
    Create a wig-format file, *outfn*, out of a pre-sorted BED or Bowtie-format
    file, *infn*.  Specify format with either filetype='bed' or
    filetype='bowtie'.
    """
    t0 = time.time()
    if outfn is None:
        fout = sys.stdout 
    else:
        fout = open(outfn,'w')

    # TODO: add some commandline options like color, trackname, visibility, etc
    # Write out the track line.  
    fout.write('track type=wiggle_0\n')
    
    for chrom, cluster_start, cluster_stop, features in clusters(infn, filetype):

        # Add 1 to cluster_start to shift WIG features so they look right in the browser (which is 1-based)
        fout.write('fixedStep chrom=%s start=%s step=1\n' % (chrom,cluster_start+1))
        cluster_len = cluster_stop - cluster_start

        # This list of zeros will be incremented at the positions of each feature
        pileup = cluster_len * [0]

        for feature_start,feature_stop in features:
            # The trick here is to subtract the cluster_start so these indices
            # go from 0 to cluster_len (instead of actual chromosomal coords),
            # since the 'fixedStep' line above indicates what chrom coord to
            # start at...
            for i in xrange(feature_start-cluster_start,feature_stop-cluster_start):
                pileup[i] += 1
        
        # Write the values (the number of stacked reads at each position) to
        # file.
        for value in pileup:
            fout.write('%s\n' % (value*options.scale))

    # Close up shop (but not if we were using stdout!)
    if outfn is not None:
        fout.close()
        
    t1 = time.time()
    sys.stderr.write('%d s elapsed\n' % (t1-t0))

if __name__ == '__main__':
    options,args = op.parse_args()
    if options.type not in dispatch_dict.keys():
        raise ValueError,'%s not a supported file type.  Choose one of %s.' % (options.type, dispatch_dict.keys())
    make_wig(options.input, options.type, options.output)
