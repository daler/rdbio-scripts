#!/usr/bin/python
usage = """
Pipeline for doing randomized intersections in parallel.  

Usage::

    ls *.bed | random-intersection-pipeline.py

Output will be a text file, randomization-results.txt, and a PDF heatmap,
heatmap.pdf, compressed into a zip file with the date and a name you provide.

**If you want to re-run, then delete the randomized-results directory.**
"""

# The meat of the randomization code is in the pybedtools package.  Currently,
# this uses a randomization that's optimized for dm3 genome.  Will take
# minor refactoring to use for others.
#
# The visualization code is in the IntersectionCluster object, separated out
# into selection, pre-sorting, clustering, labeling, and plotting.

# Built-in
import os
import sys
import optparse
import datetime

# 3rd-party
from ruffus import *
from Bio.Cluster import kcluster
from scipy import stats
from pylab import *

# my libs
import pybedtools

op = optparse.OptionParser(usage=usage)
op.add_option('--iterations',type=int, help='number of random iterations to perform')
op.add_option('--name',help='name of this run -- this is added to the final zip filename')
options,args = op.parse_args()

if not options.iterations:
    op.print_help()
    print '\nERROR: Please supply the number of iterations\n'
    sys.exit(1)

ITERATIONS = options.iterations
HEADER = ['file_a',
          'file_b',
          'count_a',
          'count_b',
          'actual',
          'lower_95th',
          'med',
          'upper_95th',
          'percentile',
          'normalized']
HEADER = '\t'.join(HEADER)+'\n'

# Get filenames from stdin
fns = []
for i in sys.stdin:
    i = i.strip()
    fns.append(i)

if len(fns) < 1:
    op.print_help()
    print '\nERROR: You need to specify some files!\n'
    sys.exit(1)

class IntersectionCluster(object):
    def __init__(self,fn,filterfunc=None, counts=False, only_antibody=False):
        """
        *fn* is a tab-delimited file containing the input data.
        """
        self.fn = fn
        self._counts = []
        self.data = csv2rec(fn, delimiter='\t')
        self._filter(filterfunc)
        self._unique_samples()
        self._names = self.samples

    def inspect(self,a,b):
        """
        Convenience function to look for *a* in self.data.file_a and *b* in self.data.file_b.
        """
        ind = []
        for i,rec in enumerate(self.data):
            if (a in rec.file_a) and (b in rec.file_b):
                ind.append(i)
        return self.data[ind]

    def _filter(self,filterfunc=None):
        """
        Resets self.data based on the return values of *filterfunc* which is
        called on each item in self.data.

        Call signature of *filterfunc* is filterfunc(rec); should return
        boolean indicating whether the rec should be kept or not.
        """
        if filterfunc is None:
            def filterfunc(rec):
                return True

        ind = []
        for i,rec in enumerate(self.data):
            if filterfunc(rec):
                ind.append(i)
       
        self.mask = array(ind)
        self.data = self.data[self.mask]
 
    def _unique_samples(self):
        """
        Converts all samples to basenames; gets the unique samples.
        """
        for i in self.data:
            i.file_a = os.path.splitext(os.path.basename(i.file_a))[0]
            i.file_b = os.path.splitext(os.path.basename(i.file_b))[0]
        self.samples = unique(self.data.file_a)

    def fill_in_matrix(self,datafunc=None,sigfunc=None):
        """
        Creates a matrix self.z which contains data from each pairwise
        intersection.  
        
        Pass an optional *datafunc* whose signature is
        datafunc(rec) if you don't want the default rec.enrichment to be the
        value in z.

        Pass an optional *sigfunc* whose signature is sigfunc(rec) and whose
        return value is True if significant or False if not.  Supply this func
        if you don't want the default of::

            def sigfunc(rec):
                if rec.percentile < 5:
                    return True
                if rec.percentile > 95:
                    return True
                return False

        """
        featurecounts = []
        if sigfunc is None:
            def sigfunc(rec):
                if rec.percentile < 0.1:
                    return True
                if rec.percentile > 99.9:
                    return True
                return False
        
        # Do the filtering etc
        self._unique_samples()

        self.z      = zeros( (len(self.samples), len(self.samples)) )
        self.sigmat = zeros( (len(self.samples), len(self.samples)) )

        # iterate through each pairwise intersection
        for ii, i in enumerate(self.samples):
            ind1 = self.data.file_a == i
            featurecounts.append(unique(self.data[ind1].count_a)[0])
            for jj, j in enumerate(self.samples):
                ind2 = self.data.file_b == j
                
                # rec1 is the record in self.data that has the data for the 
                # intersection for i vs j
                ind = ind1 & ind2                    
                rec1 = self.data[ind]

                # double check that this pairwise intersection is unique; if
                # it's not that indicates some problem with the original text
                # data file.
                assert len(rec1) == 1
                


                # the actual data we want to put in the array.
                if datafunc is None:
                    val = log2((rec1.actual+1)/(rec1.med+1))
                else:
                    val = datafunc(rec1)

                # Enter those values in z.
                row = ii 
                col = jj
                self.z[row,col] = val

                # Fill in the significance matrix
                self.sigmat[row,col] = sigfunc(rec1)

        # mask out infs and NaNs.
        self.z = ma.masked_invalid(self.z)
        self._counts = featurecounts

    def cluster(self, nclusters=None, npass=100, direction='row', initialid=None, 
                do_pca=False, add_means=True, clear_nonsig=True):
        '''Perform the clustering; results in sortinds.'''
        if nclusters is None:
            # rule of thumb for choosing k:
            # http://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set#Rule_of_thumb
            rows = self.z.shape[0]
            nclusters = int(ceil((rows/2.)**.5))

        # clear z along the diagonal
        if clear_nonsig:
            rows,cols = self.z.shape
            for r in range(rows):
                for c in range(cols):
                    if r == c:
                        self.z[r,c] = 0
            self.z[self.sigmat==0] = 0

        if do_pca:
            print 'Performing SVD on %s x %s matrix...' % (shape(self.z))
            sys.stdout.flush()
            u,s,v = svd(self.z.data)
            if direction=='row':
                u = u[:,0]
            if direction == 'col':
                u = u[0,:]
            sortind = argsort(u, kind='mergesort')
            initialid = sortind % nclusters
        
        if direction == 'row':
            transpose = 1
            axis = 1
        elif direction == 'col':
            transpose = 0
            axis = 0
        else:
            raise ValueError, 'direction must be one of "row" or "col"; %s was provided' % direction


        print 'Clustering...'
        sys.stdout.flush()
        clusterid,error,nfound = kcluster(self.z,transpose=transpose,nclusters=nclusters,npass=npass,initialid=initialid)
        
        if add_means:
            means = self.z.mean(axis=axis)
            means /= means.max()+10
            clusterid = clusterid.astype(float)
            clusterid += means

        sortind = argsort(clusterid, kind='mergesort')
        return sortind

    def plot(self, sortind, figheight=5, clear_diagonal=True, clear_nonsig=False, cmap=None, ax_rect=(.2,.2,.7,.7)):

        colorbarfrac = 0.15
        figwidth = figheight+colorbarfrac*figheight
        figsize = (figwidth,figheight)
        fig = figure(figsize=figsize)

        ax_rect = list(ax_rect)
        ax_rect[2] = ax_rect[2] - ax_rect[2]*colorbarfrac
        ax = Axes(fig,rect=ax_rect)
        fig.add_axes(ax)
        
        norm = matplotlib.colors.Normalize()
        
        # what's the normalized "0"?
        z = self.z.copy()

        norm(z)

        # Set the vmax in the colorbar to the 99th percentile to better expand
        # the colorbar
        norm.vmax = stats.scoreatpercentile(z.ravel(),99)

        zeropoint = norm(0)
        
        # split from zero to max(z) into 4 chunks.
        dcolor = (1-zeropoint)/3
        
        # various cutoffs in the colorbar I was trying...
        twofold = norm(2)
        mm = norm(abs(z.min()))
        midpoint1 = zeropoint + 1*dcolor

        # dynamically set the midpoint of the color change to the median
        # of the positive data set
        medpoint = norm(stats.scoreatpercentile(z[z>0].ravel(),50))
       

        if cmap is None:
            cdict  = {'red':  ((0.0, 0.0, 0.0),
                               (zeropoint,1.0, 1.0),
                               (medpoint, 1., 1.0),
                               (1.0, 1.0, 1.0)),

                     'green': ((0.0, 0.0, 0.0),
                               (zeropoint,1.0, 1.0),
                               (medpoint, .9, .9),
                               (1.0, 0.0, 0.0)),

                     'blue':  ((0.0, 0.0, 1.0),
                               (zeropoint,1.0, 1.0),
                               (medpoint, .0, .0),
                               (1.0, 0.0, 0.0))
                    }
            cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
            cmap.set_bad(color='k')
        
        # clear z along the diagonal
        if clear_diagonal:
            rows,cols = z.shape
            for r in range(rows):
                for c in range(cols):
                    if r == c:
                        z[r,c] = NaN
        if clear_nonsig:
            z[self.sigmat==0] = NaN
        
        z = ma.masked_invalid(z)

        pcol = ax.pcolor(z[sortind,:][:,sortind],cmap=cmap,norm=norm,edgecolors='None')
        names_sorted = self._names[sortind]
        
        nnames = len(names_sorted)
        ax.set_xticks(arange(0.5,nnames+.5))
        ax.set_xticklabels(names_sorted,rotation=90)
        ax.set_yticks(arange(0.5,nnames+.5))
        ax.set_yticklabels(names_sorted)

        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        axis('tight')
        
        cbar_padding = 0.25*colorbarfrac
        cbar_left = 1 - colorbarfrac+cbar_padding
        cbar_bottom = ax_rect[1]
        cbar_width = 0.03
        cbar_height = ax_rect[-1]
        cbar_rect = [cbar_left,cbar_bottom, cbar_width, cbar_height]
        cax = Axes(fig,rect=cbar_rect)
        fig.add_axes(cax)
        
        c = fig.colorbar(mappable=pcol,cax=cax)
        c.set_label('Enrichment score\nlog2(actual / median randomized)')
        
        counts_sorted = array(self._counts)[sortind]
        ax2 = fig.add_axes(ax.get_position(True),frameon=False)
        ax2.xaxis.tick_top()
        ax2.xaxis.set_label_position('top')
        ax2.yaxis.set_label_position('right')
        ax2.yaxis.tick_right()
        ax2.set_xticks(arange(0.5,nnames+.5))
        ax2.set_xticklabels(counts_sorted,rotation=90)
        ax2.set_yticks(arange(0.5,nnames+.5))
        ax2.set_yticklabels(counts_sorted)
        ax2.xaxis.set_ticks_position('none')
        ax2.yaxis.set_ticks_position('none')
        for t in ax2.get_xticklabels():
            t.set_rotation(90)
        ax2.axis(ax.axis())
        return fig


randomized_results_dir = 'randomized-results'
if not os.path.exists(randomized_results_dir):
    os.makedirs(randomized_results_dir)

# TASK:   Do all pairwise intersections, randomizing ITERATIONS times.
# Input:  The bed files to use
# Output: A one-line randomization output file for each pairwise comparison
def intersect_files():
    for fn1 in fns:
        for fn2 in fns:

            inputs = (fn1,fn2)
            output = '%s-%s' % (os.path.basename(fn1),os.path.basename(fn2))
            output = os.path.join(randomized_results_dir,output)
            yield (inputs,output)
@files(intersect_files)
@posttask(pybedtools.cleanup)
def random_intersections(inputs,output):
    fn1,fn2 = inputs
    print fn1,fn2
    sys.stdout.flush()
    a = pybedtools.bedtool(fn1)
    b = pybedtools.bedtool(fn2)
    results = a.randomstats(b,ITERATIONS,intersectkwargs={'u':True})
    key_order = [
        fn1,
        fn2,
        'actual',
        'lower_95th',
        'median randomized',
        'upper_95th',
        'percentile',
        'normalized',
        ]
    line = [fn1,fn2]
    line += [results[i] for i in key_order]
    line = map(str,line)
    line = '\t'.join(line)+'\n'
    fout = open(output,'w')
    fout.write(line)
    fout.close()

# TASK:   Combine all the randomized output files into one.
# Input:  All the randomization output files
# Output: A randomization report, one line for each pairwise comparison.  Also
#         has header.
def combine_randomized_files():
    inputs = []
    for _, output in intersect_files():
        inputs.append(output)
    output = 'randomization-report.txt'
    yield (inputs,output)
@files(combine_randomized_files)
@follows(random_intersections)
def combined_randomized(inputs,output):
    fout = open(output,'w')
    fout.write(HEADER)
    for f in inputs:
        fout.write(open(f).read())
    fout.close()

# TASK:   Plot a heatmap of the normalized scores.
# Input:  The randomization report
# Output: A heatmap PDF 
@follows(combined_randomized)
@files('randomization-report.txt','heatmap.pdf')
def make_heatmap(input,output):
    I = IntersectionCluster(input)
    rcParams['font.size'] = 8
    I.fill_in_matrix()
    sortind = I.cluster()
    fig = I.plot(sortind)
    fig.savefig(output)
    
# TASK:   Make a zip file of the heatmap and randomization report
# Input:  Files to zip
# Output: Zip file to send
def package_files():
    inputs = ('heatmap.pdf','randomization-report.txt')
    today = datetime.date.today()
    output = '%s-%s.zip'% (today,options.name)
    yield inputs,output
@files(package_files)
@follows(make_heatmap)
def package(inputs,output):
    cmds = ['zip',output]
    cmds.extend(inputs)
    os.system(' '.join(cmds))

pipeline_run([package],multiprocess=8)
