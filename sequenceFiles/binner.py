'''Using ideas from Jim Kent to bin features into the smallest bin they'll fit in.
1. computes start and stop arrays for bin boundaries.
2. imports a WIG file into a sqlite3 table.
3. imports a BED file into a sqlite3 table.
4. performs an intersection on them.
'''


import sqlite3
import sys
import time
import logging
import networkx as nx
from numpy import array, ceil, nonzero
import optparse
import compare_cell_types

op = optparse.OptionParser()
op.add_option('--bed',dest='bed', help='first file to import into sqlite3 db')
op.add_option('--wig',dest='wig', help='second file to import into sqlite3 db')
op.add_option('-o',dest='output',help='output file to save result as')
options,args = op.parse_args()

# Set up logging 
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('log')

starts = {}
stops = {}
start = 1
end = 5e9 # Max genome size.


divisions=4
depth=8
G = nx.balanced_tree(divisions, depth)
G2 = nx.dfs_tree(G)

def children(node):
    '''returns all childrens smaller than node.'''
    L = []
    for i in nx.dfs_successor(G2,node).values():
        L.extend(i)
    return L

# generate labels for tree.
nodes_per_level = [divisions**i for i in range(depth+1)]
nnodes = sum(nodes_per_level)    
node = 0

#step_sizes = [512e6, 64e6, 8e6, 1e6, 128e3]
#nodes_per_level = [int(ceil(end/i)) for i in step_sizes]


for i in nodes_per_level:
    step = int(ceil(end/float(i)))
    start = 1
    for k in range(i):
        stop = start + step
        starts[node] = start
        stops[node] = stop
        node += 1
        start += step

stopdata = array(stops.values())
startdata = array(starts.values())


"""
# bin sizes empirically determined from Genome Browser paper
bin_sizes = [end, 512e6, 64e6, 8e6, 1e6, 128e3, 16e3]

# based on bin sizes, this is the nuber of bins in each "level"
bins_per_level = [int(ceil(end/i)) for i in bin_sizes]

bin = 0
for i in bins_per_level:
    step = int(ceil(end/float(i)))
    start = 1
    for k in range(i):
        stop = start + step
        starts[bin] = start
        stops[bin] = stop
        bin += 1
        start += step
stopdata = array(stops.values())
startdata = array(starts.values())
"""
def determine_bin(start,stop):
    '''Returns the bin number of the smallest bin the feature
    will fit within and a list of bins that encompass this smallest bin.'''
    inds = nonzero( (start >= startdata) & (stop <= stopdata))[0]
    try:
        ind = inds[-1]
    except IndexError:
        raise ValueError, 'No bins will contain this feature!'
    return ind, inds.tolist()

# initialize db
conn = sqlite3.connect('test.db')
c = conn.cursor()
c.execute('drop table if exists features1')
c.execute('drop table if exists features2')
c.execute('create table if not exists features1 (start INTEGER, stop INTEGER, bin INTEGER, value FLOAT, chrom TEXT, span INTEGER)')
c.execute('create table if not exists features2 (start INTEGER, stop INTEGER, bin INTEGER, chrom TEXT)')

logger.info('importing WIG file into database...')
#wigs = compare_cell_types.wig_features('tests/compare_cell_types/wig.wig')
wigs = compare_cell_types.wig_features(options.wig)
for w in wigs:
    bin = determine_bin(w.start, w.stop)[0]
    c.execute('insert into features1 (chrom,start,stop,value,bin,span) VALUES (?,?,?,?,?,?)',
              (w.chrom, w.start, w.stop, w.value, bin, w.span))
conn.commit()

#bedfeatures = compare_cell_types.bed_features('tests/compare_cell_types/bed.bed')
bedfeatures = compare_cell_types.bed_features(options.bed)

# assume that bedfeatures is sorted by chromosome!
wigfeatures = []
for b in bedfeatures:
    bin, larger_bins = determine_bin(b.start,b.stop)
    
    #all_bins = larger_bins + children(bin)
    all_bins = larger_bins

    logger.info('%s,\n\t%s, %s' % (b,bin,all_bins))

    # in_str will be something like "?1,?2,?3,?4"
    in_str = ['?%s,' % i for i in range(1,len(all_bins)+1)]
    in_str = ''.join(in_str)
    in_str = in_str[:-1]
    data = all_bins
    data.append(b.chrom)
    candidates = c.execute('select chrom,start,stop,bin,value,span from features1 where bin in (%s) and chrom=?' % in_str,
                            data).fetchall()    

    logger.info('found %s possible WIG features' % len(candidates))
    starts = array([i[1] for i in candidates])
    stops = array([i[2] for i in candidates])
    try:
        inds = nonzero((starts >= b.start) & (stops <= b.stop))[0]
    except IndexError:
        inds = [nonzero((starts >= b.start) & (stops <= b.stop))]
    logger.info('\tof these, %s overlap' % len(inds))
    for ind in inds:
        chrom,start,stop,bin,value,span = candidates[ind]
        f = compare_cell_types.feature(start=start,chrom=chrom,stop=stop,value=value,span=span)
        wigfeatures.append(f)
#compare_cell_types.write_wig('overlapping-wigs.wig',wigfeatures)
compare_cell_types.write_wig(options.output,wigfeatures)
conn.close()


#888888888888888888888888888888888888888888888
# test functions
#888888888888888888888888888888888888888888888
def test_binning():
    '''Not really a test but more of a debugging tool.'''
    featurestart = int(1e6-100)
    featurestop =  int(1e6-90)
    inds = nonzero((featurestart>=startdata) & (featurestop<=stopdata))[0]
    ind = inds[-1]
    print 'feature is at   %s - %s' % (featurestart,featurestop)
    print 'smallest bin is',startdata[ind],'-',stopdata[ind]
    print 'bins to search:'
    for i in inds:
        print '%s (%s,%s)' % (i,startdata[i],stopdata[i])

def test_binning2():
    '''Generates 100 random start and stop coordinates and runs the determine_bin
    function on them.  Asserts that '''
    import random
    for k in range(100):
        start = random.randint(1,end)
        stop = random.randint(1,end)
        ind, inds = determine_bin(start,stop)
        assert ind is not None
        assert stop < stopdata[ind]
        assert start > startdata[ind]
        for i in inds:
            assert stop < stopdata[i]
            assert start > startdata[i]

if __name__ == "__main__":
    pass
