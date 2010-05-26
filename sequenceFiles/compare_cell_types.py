#!/usr/bin/python

__doc__ = '''
this script takes an input WIG file and a corresponding BED file.
For each feature in the BED file, it extracts the WIG data and sends 
the output to a new file.
'''
import gzip
import os
import sys


class feature(object):
    def __init__(self,chrom,start,stop,value=None,span=None):
        '''
        Simple feature class.

        chrom : str
            Chromosome name

        start : int
            Integer start coordinate

        stop : int
            Integer stop coordinate

            '''
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.value = value
        self.span = span

    def __repr__(self):
        return 'Feature object, %s:%s-%s' % (self.chrom, self.start, self.stop)

    def within(self,other,equal=True):
        '''Test to see whether this feature is within the provided range.
        '''
        if self.chrom != other.chrom:
            return False
        if self.start > other.start:
            if self.stop < other.stop:
                return True
        return False

def bed_features(fn,first=None,last=None,equal=False):
    '''Generator function that returns features from BED file.'''
    if os.path.splitext(fn)[1] == '.gz':
        f = gzip.open(fn)
    else:
        f = open(fn)
    for line in f:
        if 'track' in line:
            continue
        L = line.rstrip().split('\t')
        if len(L)<=1:
            continue
        start = int(L[1])
        stop = int(L[2])
        if first is not None and start < first:
            continue
        if last is not None and stop > last:
            return
        yield feature(chrom=L[0], start=start, stop=stop)

def wig_features(fn,first=None,last=None,equal=False):
    '''Generator function that returns features from WIG file.'''
    
    if os.path.splitext(fn)[1] == '.gz':
        f = gzip.open(fn)
    else:
        f = open(fn)
    track_data = {}
    for line in f:
        L = line.rstrip().split()
        try: # determine if this is track description or data
            int(L[0])
            valid_data = True
        except ValueError:
            valid_data = False
        
        if valid_data:
            try:
                track_data
            except NameError:
                print 'No header found.  Exiting.'
                sys.exit()
            start = int(L[0])
            value = float(L[1])
            chrom = track_data['chrom']
            stop = start + track_data['span']
            if first is not None and start <= first:
                continue
            if last is not None and stop >= last:
                return
            yield feature(chrom=chrom, start=start, stop=stop, value=value, span=track_data['span'])
            
        elif not valid_data:
            for i in L:
                if '=' in i:
                    key, value = i.split('=')
                    try:
                        track_data[key] = int(value)
                    except ValueError:
                        track_data[key] = value

def write_wig(output_fn, wigs):
    '''If output_fn is an open file, then use that.  Otherwise, open a new file for writing.'''
    if type(output_fn) is str:
        fout = open(output_fn,'w')
    else:
        fout = output_fn
    chrom = None
    for w in wigs:
        if w.chrom != chrom:
            chrom = w.chrom
            fout.write('type=wiggle_0\n')
            fout.write('variableStep chrom=%s span=%s\n' % (w.chrom,w.span))
        fout.write('%s\t%s\n' % (w.start,w.value))

    # if input was a file handle, leave it open
    if type(output_fn) is str:
        fout.close()

            
if __name__ == "__main__":   

    import sys
    import optparse
    op = optparse.OptionParser()
    op.add_option('--bed', dest='bedfn', help='input bed file')
    op.add_option('--wig', dest='wigfn', help='input wig file')
    op.add_option('--output', dest='outputfn', help='output fn, will be a wig file')
    options,args = op.parse_args()

    if len(sys.argv)>1:
        write_wig(options.outputfn, intersect(options.bedfn,options.wigfn))
    else:
        op.print_help()
