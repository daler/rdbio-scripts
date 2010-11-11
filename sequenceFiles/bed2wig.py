#!/usr/bin/python

import optparse
import time
import sys

"""
Script to convert SAM, BED or Bowtie-formatted files into a "pileup" or
"coverage" wig file.

***IMPORTANT***

    Assumes that the input file is sorted by chromosome, then by start
    position.  Run this script with the -h option for info on how to sort.

How does it work?
-----------------
The strategy is to bin together contiguously overlapping features into clusters
of reads, computing the cumulative coverage on each cluster.  Then the cluster
is written to file as WIG format and removed from memory. The result is that
only one cluster is in memory at a time.  

In contrast to previous implementations that intitialized an array the size of
the genome, long stretches of empty genome are simply not checked here -- only
the features in the input file are used.

Other formats
-------------
Support for other input file types can be implemented by writing another
iterator that returns chrom, start, stop (start and stop must be integers) and
adding that to the dispatch dictionary.  

Revision history
----------------
Created Jan 2010 Ryan Dale

Modified Mar 8 2010 Ryan Dale:
    * added support for SAM format

Modified Apr 15 2010 Ryan Dale:
    * added support for scale factor

Modified Aug 2 2010 Ryan Dale:
    * "alwaysZero" added to output wig to show the zero values when autoscaling
    * added help text for SAM and BAM file sorting
    * added example usage
    * sys.exit with no input args
    * added option checking with meaningful error messages

Modified Aug 9 2010 Ryan Dale:
    * added --strand option; changed iterators to return strand as well
    * added --track option to add titles, etc

Modified Aug 19 2010 Ryan Dale:
    * added piping support -- make_wig, cluster, and all iterators use open file handles now
"""

usage = """
Convert a pre-sorted SAM, BED or Bowtie file into WIG format.

***IMPORTANT***

    Assumes that the input file is sorted by chromosome, then by start
    position.  
    
    For a *Bowtie* format file, sort your file with:

        sort -k 3,3 -k 4n INPUTFILE > OUTPUTFILE

    For a *BED* file, use::

        sort -k 1,1 -k 2n INPUTFILE > OUTPUTFILE

    For a *SAM* file, use something like the following (this example uses
    samtools to first filter out unmappable reads)::
        
        samtools view -S -F 0x0004 INPUTFILE | sort -k 3,3 -k 4n > OUTPUTFILE

    If all you have is a *BAM* file, then you can sort with samtools first,
    then output as SAM:

        samtools sort INPUTFILE SORTED
        samtools view -F 0x0004 SORTED.bam > OUTPUTFILE


Example usage:
--------------

Minimal usage:
    
    python bed2wig.py -i SORTED_SAM_FILE --type sam > out.wig

Using other options: 

    python bed2wig.py -i in.bed --type bed -o out.wig --scale 1.5 --strand + --track 'name="Control tag density" color=128,0,0'

"""

op = optparse.OptionParser(usage=usage)
op.add_option('-i',     dest='input',  help='Required input file to be converted to WIG format. '
                                            'Can be BED, SAM, or Bowtie format')
op.add_option('--type', dest='type',   help='Required type of input file, one of "bed", "sam" or "bowtie".')
op.add_option('-o',     dest='output', help='Optional output file to write to. If unspecified, writes to stdout.')
op.add_option('--scale',dest='scale',  help='Optional scale factor to multiply all WIG values by. Useful '
                                            'for comparing datasets with differing library sizes', type=float,
                                            default=1.0)
op.add_option('--strand',dest='strand',help='Strand to use for output.  Can be "+" or "-".'
                                            'Default is to ignore strand.',
                                            default='.')
op.add_option('--track',dest='track',help='Additional track info to write, e.g, \'name="Input, +" color=128,0,0\'. '
                                          'Be sure to use nested or escaped quotes if your track names have spaces in them!',
                                     default='')
op.add_option('--verbose',action='store_true',help='Print progress to stderr')

# NOTE: Support for other input file types can be implemented by writing
# another iterator that returns chrom, start, stop (start and stop must be
# integers) and adding that to the dispatch dictionary.

def samfile_iterator(fn):
    """
    Iterates through a SAM file *fn*, returning chrom,start,stop for each
    line where there was a mapped read.  start and stop are converted to ints.

    Doesn't do anything fancy with CIGAR strings -- just adds len(seq) to start
    to get the stop coord.

    Ignores unmapped reads.
    """
    for line in fn:
        if line.startswith('@'):
            continue
        L = line.strip().split('\t')
        flag = int(L[1])
        if flag & 0x0004:
            continue
        strand = '+'
        if flag & 0x0010:
            strand = '-'
        chrom = L[2]
        start = int(L[3])-1
        stop = start + len(L[9])
        yield (chrom,start,stop,strand)
            

def bedfile_iterator(fn):
    """
    Iterates through a BED file *fn*, returning chrom,start,stop for each line
    (start and stop are converted to int)

    Track lines are gracefully ignored.
    """
    for line in fn:
        L = line.strip().split()
        chrom = L[0]
        try:
            start = int(L[1])
        except:
            if line.startswith('track'):
                continue
            if line.startswith('browser'):
                continue
            else:
                raise ValueError,'mis-formatted bed line:\n%s'%line
        stop = int(L[2])
        try:
            strand = L[5]
        except IndexError:
            strand = '.'
        yield (chrom,start,stop,strand)

def bowtie_iterator(fn):
    """
    Iterates through a Bowtie-format file *fn*, returning chrom,start,stop for
    each line (start and stop are converted to int)
    """
    for line in fn:
        L = line.strip().split()
        chrom = L[2]
        strand = L[1]
        start = int(L[3])
        stop = start + len(L[4])-1
        yield (chrom,start,stop,strand)

# Makes the connection between filetype name and the iterator to use for that
# filetype
dispatch_dict = {
                 'bed':bedfile_iterator,
                 'bowtie':bowtie_iterator,
                 'sam':samfile_iterator,
                }

def clusters(infn, filetype, use_strand='.',verbose=False):
    """
    Yields clusters of overlapping reads along with the chromsome and cluster boundaries.

    *strand* is one of '+','-' or '.'

    Return value is of the form (chrom, cluster_start, cluster_stop, features)
    where *features* is a list of (start,stop) tuples.
    """

    # Decide which iterator to use for this filetype.
    try:
        iterator = dispatch_dict[filetype](infn)
    except KeyError:
        raise ValueError, "Input filetype %s not supported; only %s currently supported" % (filetype,dispatch_dict.keys())

    # If there's nothing in the iterator, then don't return the last cluster
    # (since there won't be one)
    ret_last = False
    last_chrom = None
    for chrom,start,stop,strand in iterator:
        # If we enter this loop, there's something in the iterator and so we
        # will eventually have to return the last cluster. 
        if use_strand != '.':
            if strand != use_strand:
                continue

        ret_last = True

        # If we're on a new chromosome, reset everything.
        if chrom != last_chrom:
            if verbose:
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
            

    if ret_last:
        # Return the last cluster (since it won't be yielded by the else-clause above)
        yield chrom, cluster_start, cluster_stop, features

def make_wig(infn, filetype, outfn=None, use_strand='.',trackinfo=None,verbose=False):
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
    if trackinfo is None:
        trackinfo = ''
    # TODO: add some commandline options like color, trackname, visibility, etc
    # Write out the track line.  
    fout.write('track type=wiggle_0 alwaysZero=on %s\n' % trackinfo)

    
    for chrom, cluster_start, cluster_stop, features in clusters(infn, filetype, use_strand, verbose):
        if chrom is None:
            continue
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
    if verbose:
        sys.stderr.write('%d s elapsed\n' % (t1-t0))

if __name__ == '__main__':
    options,args = op.parse_args()
    def operr(msg):
        op.print_help()
        sys.stderr.write('\n***ERROR: %s***\n' % msg)
        sys.exit(1)

    if not options.type:
        operr('Please specify a type!')
    if options.type not in dispatch_dict.keys():
        operr('"%s" is not a supported file type.  Choose one of %s.\n' % (options.type, dispatch_dict.keys()))
    if not options.input:
        operr('Please specify an input file!')

    if options.input == 'stdin':
        input_handle = sys.stdin
    else:
        input_handle = open(options.input)
    make_wig(input_handle, options.type, options.output, options.strand, options.track,options.verbose)
