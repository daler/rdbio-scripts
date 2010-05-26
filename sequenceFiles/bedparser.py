"""
Module with parsers for BED, WIG, BEDGRAPH, and SAM format files.
"""
import gzip
import pdb
import os
import shlex

def parsetrackline(trackline):
    """Parses a trackline into key/value pairs which are converted into a
    dictionary."""
    assert trackline.startswith('track ')
    t = trackline[6:]
    items = shlex.split(t) # shlex preserves quoted substrings
    keys = []
    values = []
    for i in items:
        key,value = i.split('=')
        if " " in value:
            value = '"%s"' % value
        keys.append(key)
        values.append(value)
    return dict(zip(keys,values))

class Track(object):
    """Track lines with attributes as defined here:
    
        http://genome.ucsc.edu/goldenPath/help/customTrack.html

    Converts anything that can be into a float.  Otherwise a string.  All defaults
    are None.
    """
    def __init__(self, **kwargs):
        self._validnames = ['name',
                            'description',
                            'visibility',
                            'color',
                            'useScore',
                            'priority',
                            'offset',
                            'url',
                            'db',
                            'itemRgb',
                            'group']

        for attr in self._validnames:
            try:
                value = kwargs[attr]
                
                if attr not in  ['name', 'description', 'group','db']:
                    # don't want to convert to float a track name.
                    try:
                        value = float(value)
                    except ValueError:
                        pass
            except KeyError:
                value = None
            setattr(self,attr,value)

        for key,value in kwargs.iteritems():
            if key not in self._validnames:
                raise Warning('%s is not a valid track line variable' % key)
            
    def __repr__(self):
        hasanyvalue = False
        s = 'track '
        for attr in self._validnames:
            value = getattr(self,attr)
            if value is not None:
                hasanyvalue = True
                s += attr+'='+str(value)+' '
        if not hasanyvalue:
            return 'track name="User track"'
        else:
            return s

class bedfeature(object):
    def __init__(self, chr,start,stop,
                 name=None,value=None,strand=None,
                 thickStart=None,thickStop=None,itemRGB=None,
                 blockCount=None,blockSizes=None,blockStarts=None, 
                 track=Track()):
        self.chr=chr
        self.start=int(start)
        self.stop=int(stop)
        self.name=name
        self.strand=strand
        self.track=track
        try:
            self.value=float(value)
        except TypeError:
            self.value=value
        try:
            self.thickStart=int(thickStart)
        except TypeError:
            self.thickStart=thickStart
        try:
            self.thickStop=int(thickStop)
        except TypeError:
            self.thickStop=thickStop
        try:
            self.blockCount=int(blockCount)
        except TypeError:
            self.blockCount=blockCount
        
        self.itemRGB=itemRGB
        self.blockSizes=blockSizes
        self.blockStarts=blockStarts

    def __repr__(self):
        return 'bed feature: %s:%s-%s' % (self.chr,self.start,self.stop)
    
    def tostring(self):
        """Prints the bed record suitable for writing to file, newline included.
        
        In the interest of speed, does not do error-checking.
        """
        items = [self.chr, 
                 self.start, 
                 self.stop, 
                 self.name, 
                 self.value, 
                 self.strand, 
                 self.thickStart,
                 self.thickStop, 
                 self.itemRGB,
                 self.blockCount,
                 self.blockSizes, 
                 self.blockStarts]
        printables = []
        for item in items:
            if item is None:
                printables.append('')
            else:
                printables.append(str(item))
        
        return '\t'.join(printables).rstrip()+'\n'

class wigfeature(object):
    def __init__(self, chr ,start,value,span):
        self.chr=chr
        self.start=int(start)
        self.value=float(value)
        self.stop=self.start + int(span)
    def __repr__(self):      
        return 'wig feature: %s:%s-%s' % (self.chr,self.start,self.stop)

class bedgraphfeature(object):
    def __init__(self, chr,start,stop,value):
        self.chr=chr
        self.start=int(start)
        self.stop=int(stop)
        self.value=float(value)
    def __repr__(self):
        return 'bedgraph feature: %s:%s-%s %s' % (self.chr,self.start,self.stop, self.value)


class bedfile(object):
    """Iterator object, with __iter__ defined, that moves through
    features in a BED-format file.  A new bedfeature object is 
    created for each line.
    
    Usage::

        for feature in bedfile('a.bed'):
            print feature.chr
            print feature.start
            print feature.stop
            print 'length:', feature.stop-feature.start
        """
    def __init__(self,f):
        if type(f) is str:
            self.stringfn = True
            if os.path.splitext(f)[-1] == '.gz':
                self.file = gzip.open(self.fn)
            else:
                self.file = open(f)
        else:
            self.stringfn = False
            self.file = f

    def __iter__(self):
        f = self.file
        track = Track()
        for line in f:
            L = line.rstrip().split('\t')
            if 'track' in line:
                tracklinedict = parsetrackline(line)
                track = Track(**tracklinedict)
                continue

            if 'browser' in line:
                continue

            # blank line?
            if len(line.rstrip()) == 0:
                continue
            #assert len(L) >= 3
            args = [None for i in range(12)]
            args[:len(L)] = L
            args.append(track)
            yield bedfeature(*args)
        if self.stringfn:
            f.close()

    def __repr__(self):
        return 'bedfile object with %s features. file=%s' % (self.count, self.fn)

class wigfile(object):
    
    def __init__(self,fn):
        self.fn = fn
    
    def __iter__(self):
        if os.path.splitext(self.fn)[-1] == '.gz':
            f = gzip.open(self.fn)
        else:
            f = open(self.fn)
        for line in f:
            if line.startswith('#'):
                continue
            if 'track' in line:
                continue

            # next chromosome, so parse this line.
            if 'chrom=' in line:
                assert 'variableStep' in line 
                # (not sure how to deal with other types besides
                # variableStep yet)

                L = line.rstrip().split(' ')
                for i in L:
                    if 'chrom' in i:
                        chrom=i.split('=')[1]
                    if 'span' in i:
                        span=i.split('=')[1]
                continue
            
            # made it this far?  You're on a data line.
            L = line.rstrip().split('\t')
            assert len(L) == 2
            start = L[0]
            value = L[1]
            yield wigfeature(chrom,start,value,span)
        f.close()

        
    def __repr__(self):
        return 'wigfile object, file=%s' % self.fn

class bedgraph(object):
    
    def __init__(self,fn):
        self.fn = fn
    
    def __iter__(self):
        if os.path.splitext(self.fn)[-1] == '.gz':
            f = gzip.open(self.fn)
        else:
            f = open(self.fn)
        for line in f:
            if line.startswith('#'):
                continue
            if 'track' in line:
                assert 'bedGraph' in line
                continue

            L = line.rstrip().split('\t')
            assert len(L) == 4
            chr = L[0]
            start = L[1]
            end = L[2]
            value = L[3]
            yield bedgraphfeature(chr,start,end,value)
        f.close()
    

    def __repr__(self):
        return 'bedgraph object, file=%s' % self.fn 

# extremely naive!
class samfeature(object):
    def __init__(self, chr, start, stop, strand):
        self.chr=chr
        self.start=int(start)
        self.stop=int(stop)
        self.strand=strand
    def __repr__(self):
        return 'SAM feature: %s:%s-%s (%s)' % (self.chr,self.start,self.stop,self.strand)
# extremely naive.
class samfile(object):
    def __init__(self,fn):
        self.fn = fn
    def __iter__(self):
        if os.path.splitext(self.fn)[-1] == '.gz':
            f = gzip.open(self.fn)
        else:
            f = open(self.fn)
        for line in f:
            L = line.rstrip().split('\t')
            name,flag,rname,pos,mapq,cigar,mrnm,mpos,isize,seq,qual = L[0:11]
            chr = rname
            flag = int(flag)
            start = int(pos)
            stop = start + len(seq)
            if flag == 16:
                strand = '-'
            if flag == 0:
                strand = '+'
            else: 
                #raise ValueError, 'flag is %s' % flag
                #pdb.set_trace()
                pass
            yield samfeature(chr,start,stop,strand)
        f.close()
    def __repr__(self):
        return 'bedfile object with %s features. file=%s' % (self.count, self.fn)

class gfffile(object):
    """Iterator object, with __iter__ defined, that moves through
    features in a GFF-format file.  A new gfffeature object is 
    created for each line.
    
    Usage::

        for feature in gfffile('a.bed'):
            print feature.chr
            print feature.start
            print feature.stop
            print feature.featuretype
            print feature.desc
            print 'length:', feature.stop-feature.start
        """
    def __init__(self,f,strvals=False):
        if type(f) is str:
            self.stringfn = True
            if os.path.splitext(f)[-1] == '.gz':
                self.file = gzip.open(self.fn)
            else:
                self.file = open(f)
        else:
            self.stringfn = False
            self.file = f
        self.strvals = strvals

    def __iter__(self):
        f = self.file
        for line in f:
            line = line.rstrip()
            if line.startswith('#') or len(line) == 0:
                continue
            L = line.rstrip().split('\t')
            args = [None for i in range(9)]
            args[:len(L)] = L
            args.append(self.strvals)
            yield gfffeature(*args)
        if self.stringfn:
            f.close()

    def __repr__(self):
        return 'gfffile object (file=%s)' % (self.file)

class gfffeature(object):
    
    class attributes(object):
        def __init__(self):
            self._attrs = []  # will hold a list of attributes added to the object.

    def __init__(self, chr, source, featuretype, start, stop,
                 value,strand,phase,attributes,strvals=False):
        
        if not strvals:
            self.chr=chr
            self.source=source
            self.featuretype=featuretype
            self.start=int(start)
            self.stop=int(stop)
            try:
                self.value=float(value)
            except ValueError:
                self.value=None
            self.strand=strand
            try: 
                self.phase=int(phase)
            except ValueError:
                self.phase=None
        if strvals: 
            self.chr=chr
            self.source=source
            self.featuretype=featuretype
            self.start=start
            self.stop=stop
            self.value=value
            self.strand=strand
            self.phase=phase
            
        self._strattributes = attributes # save these for later printing out.
        # parse description
        self.attributes = gfffeature.attributes()
        items = attributes.split(';')
        for item in items:
            if len(item) > 0:
                field,value = item.split('=')
            values = value.split(',')
            setattr(self.attributes,field,values)
            self.attributes._attrs.append(field)

    def __repr__(self):
        return 'GFF %s feature: %s:%s-%s' % (self.featuretype,self.chr,self.start,self.stop)
    
    def tostring(self):
        """Prints the GFF record suitable for writing to file, newline included.
        
        In the interest of speed, does not do error-checking.
        """
        attributes = ''
        for attr in self.attributes._attrs:
            values = getattr(self.attributes,attr)
            attributes += attr+'='+','.join(values)+';'

        items = [self.chr, 
                 self.source,
                 self.featuretype,
                 self.start, 
                 self.stop, 
                 self.value, 
                 self.strand, 
                 self.phase,
                 attributes]
        printables = []
        for item in items:
            if item is None:
                printables.append('.')
            else:
                printables.append(str(item))
        return '\t'.join(printables).rstrip()+'\n'

if __name__ == "__main__":
    import sys
    for i in gfffile(sys.argv[1],strvals=True):
        pass
