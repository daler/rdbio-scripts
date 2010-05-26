#!/usr/bin/python
"""



Usage example::

"""
class Cluster(object):
    
    def __init__(self,feature=None, forceclustersize=None, scorefunc=None, minclustersize=None, minclusterscore=None, minfeaturecount=None, gapwidth=None, threshold=None):
        
        self.forceclustersize = forceclustersize
        self.minclustersize = minclustersize
        self.minfeaturecount = minfeaturecount
        self.gapwidth = gapwidth
        self.threshold = threshold
        self.minclusterscore = minclusterscore
        self.clusterscores = [] 
        self.value = 0
        if scorefunc is None:
            self.scorefunc = sum

        if feature is not None:
            self.chrom = feature.chrom
            self.start = feature.start
            self.stop = feature.stop
            self.values = [feature.value]
            self.features = [feature]

        else:
            self.chrom=None
            self.start=None
            self.stop=None
            self.values = []
            self.features = []

    @property
    def count(self):
        return len(self.features)
   
    @property
    def clusterscore(self):
        return self.scorefunc(self.clusterscores)

    def check_yield(self):
        """
        Checks to see if this cluster should be returned.
        """
        if self.minclustersize and (len(self) < self.minclustersize):
            return False
        
        if self.minfeaturecount and (self.count < self.minfeaturecount):
            return False

        if self.forceclustersize and (len(self) < self.forceclustersize):
            self.stop = self.start + self.forceclustersize - 1

        if self.minclusterscore and (self.clusterscore < self.minclusterscore):
            return False

        return True 
            

    def check_feature(self,feature):
        """
        Checks to see if adding this *feature* will still make this
        a valid cluster.

        Adding a feature could invalidate the cluster by:

          * making the cluster too big (go over forceclustersize)
          * not passing the threshold
          * creating too big of an internal gap
          * being on a different chrom
        """
        
        # If there's nothing in this cluster yet, then adding the feature won't invalidate it
        if (self.start is None) and (self.stop is None) and (self.chrom is None):
            return True

        # different chrom would make it an invalid cluster
        if (self.chrom) and (feature.chrom != self.chrom):
            return False

        # under threshold would make it an invalid cluster
        if feature.value < self.threshold:
            return False

        # too big of a gap makes it invalid
        if feature.start - self.stop >= self.gapwidth:
            return False
     
        # enlarging it past forceclustersize makes it invalid
        if self.forceclustersize:
            if feature.stop - self.start > self.forceclustersize:
                return False

        return True

    def unique_features(self):
        """Returns the number of unique features, defined as the length of the
        set of all by chr,start,stop,strand tuples."""
        tuples = []
        for f in self.features:
            tuples.append( (f.chrom,f.start,f.stop,f.strand) )
        return len(set(tuples))

    def extend(self,feature):
        """
        Does not check for validity -- you should call self.check_feature() to test
        for this.
        """
        if not self.start:
            self.start = feature.start
        self.stop = feature.stop
        self.features.append(feature)
        self.clusterscores.append(feature.value)

    def __len__(self):
        if self.start:
            return self.stop-self.start
        else:
            return 0
    
    def tostring(self):
        return '%s\t%s\t%s\t%s\t%s\t%s\n'%(self.chrom,self.start,self.stop,'.',self.value,'+')
    
class Feature(object):

    def __init__(self,chrom,start,stop,value,strand='+'):
        self.chrom=chrom
        self.start=start
        self.stop=stop
        self.value=value
        self.strand=strand

    def tostring(self):
        return '%s\t%s\t%s\t%s\t%s\t%s\n'%(self.chrom,self.start,self.stop,'.',self.value,self.strand)

def bed_iterator(fn,forcevalue=None):
    
    for line in open(fn):
        if 'track' in line:
            continue
        if 'browser' in line:
            continue
        L = line.strip().split()
        chrom,start,stop = L[:3]
        if not forcevalue:
            try:
                v = float(L[4])
            except ValueError:
                raise ValueError, 'Need to specify a value since there is not one in the bed file.'  
        else:
            v = forcevalue
        yield Feature(chrom, int(start), int(stop),v)


def brennecke_cluster(fn):
    kwargs = {'gapwidth':  1e15,
              'threshold': 0,
              'minclusterscore':5,
              'forceclustersize': 5000}

    fn1 = fn+'.bren.clusters'
    fout = open(fn1,'w')
    fout.write('track name="5-kb windows"\n')
    cluster = Cluster(**kwargs)
    for feature in bed_iterator(fn):
        if cluster.check_feature(feature):
            cluster.extend(feature)
        else:
            if cluster.check_yield():
                if cluster.chrom:
                    cluster.value = cluster.clusterscore
                    fout.write(cluster.tostring())
            cluster = Cluster(feature=feature, **kwargs)
    fout.close() 
 
    kwargs = {'gapwidth':  20000,
              'threshold': 0,
              'minfeaturecount': None,
              'minclustersize': None,
              'forceclustersize': None}
    
    fn2 = fn1+'.clustered'
    fout = open(fn2,'w')
    fout.write('track name="clustered"\n')
    cluster = Cluster(**kwargs)
    for feature in bed_iterator(fn1):
        if cluster.check_feature(feature):
            cluster.extend(feature)
        else:
            if cluster.check_yield():
                if cluster.chrom:
                    fout.write(cluster.tostring())
            cluster = Cluster(feature=feature, **kwargs)
    fout.close()
   

def hannon_cluster(fn):
    unique_features = 3

    kwargs = {'gapwidth':  1e15,
              'threshold': 1,
              'minfeaturecount': 3,
              'minclustersize': 1,
              'forceclustersize': 200}
    fn1 = fn+'.clusters'
    fout = open(fn1,'w')
    cluster = Cluster(**kwargs)
    for feature in bed_iterator(fn,2):
        if cluster.check_feature(feature):
            cluster.extend(feature)
        else:
            if cluster.check_yield():
                if cluster.unique_features() >= unique_features:
                    cluster.value = cluster.unique_features()
                    fout.write(cluster.tostring())
            cluster = Cluster(feature=feature, **kwargs)
    fout.close() 

    kwargs = {'gapwidth':  200,
              'threshold': 1,
              'minfeaturecount': None,
              'minclustersize': None,
              'forceclustersize': None}
    
    fn2 = fn1+'.clustered'
    fout = open(fn2,'w')
    cluster = Cluster(**kwargs)
    for feature in bed_iterator(fn1):
        if cluster.check_feature(feature):
            cluster.extend(feature)
        else:
            if cluster.check_yield():
                sumvalues = sum([f.value for f in cluster.features])
                if sumvalues > 40:
                    density = float(sumvalues) / len(cluster)
                    cluster.clusterscores = [density]
                    fout.write(cluster.tostring())
            cluster = Cluster(feature=feature, **kwargs)
    fout.close()

if __name__ == "__main__":
    import sys
    fn = sys.argv[1]
    
    try:
        other = sys.argv[2]
    except IndexError:
        other = None

    if other == 'hannon':
        hannon_cluster(fn)

    if other == 'brenn':
        brennecke_cluster(fn)
