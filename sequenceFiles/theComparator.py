#!/usr/env/python

"""
This script does various comparisons between two input BED files.

Requires BEDTools to be available on the path.
"""

import os, optparse, sys
import tempfile
from pylab import *

op = optparse.OptionParser()
op.add_option('-a', dest='a', help='First BED file to compare')
op.add_option('-b', dest='b', help='Second BED file to compare')
op.add_option('-w', dest='w', default='0', help='Slop to use when determining overlap')
options,args = op.parse_args()

def do_overlap(a,b):
    lines_a = 0
    for line in open(a):
        lines_a += 1
    
    lines_b = 0
    for line in open(b):
        lines_b += 1
    label_a = os.path.basename(a)
    label_a = label_a.split('_')[0]
    label_b = os.path.basename(b)
    label_b = label_b.split('_')[0]
    overlaps_fn = 'features-in-%s-overlapping-%s-window-%s.bed' % (label_a,label_b,options.w)
    trackname = '"overlap between %s and %s"' % (label_a,label_b)
    cmds = ['windowBed',
            '-a', a,
            '-b', b,
            '-w', options.w,
            '-u',
            '|' 'bedTracknamer.py',
            '-i','stdin',
            '--trackname', trackname,
            '>', 
            overlaps_fn]
    os.system(' '.join(cmds))
    
    # how many overlaps were there?
    N_overlaps = 0
    for line in open(overlaps_fn):
        N_overlaps += 1
    print '\n%s out of %s features in %s overlapped features in %s' % (N_overlaps,
                                                                     lines_a,
                                                                     label_a,
                                                                     label_b)

    non_overlaps_fn = 'features-in-%s-that-do-not-overlap-%s-window-%s.bed' % (label_a,label_b,options.w)
    trackname = '"in %s but not %s"' % (label_a, label_b)
    cmds = ['windowBed',
            '-a', a,
            '-b', b,
            '-w', options.w,
            '-v',
            '|','bedTracknamer.py',
            '-i','stdin',
            '--trackname',trackname,
            '>', 
            non_overlaps_fn]
    os.system(' '.join(cmds))
     
    # how many overlaps were there?
    N_overlaps = 0
    for line in open(non_overlaps_fn):
        N_overlaps += 1
    print '\n%s out of %s features in %s did not overlap any features in %s' % (N_overlaps,
                                                                     lines_a,
                                                                     label_a,
                                                                     label_b)


def do_proximity(a,b):

    proximity_fn = tempfile.mktemp()
    cmds = ['closestBed',
            '-a',a,
            '-b',b,
            '-t','"first"',
            '>',proximity_fn]
    os.system(' '.join(cmds))
    proximities = []
    for line in open(proximity_fn):
        values = line.strip().split('\t')
        chra = values[0]
        starta = int(values[1])
        stopa = int(values[2])
        chrb = values[3]
        startb = int(values[4])
        stopb = int(values[5])
        centera = starta + ((stopa-starta)/2)
        centerb = startb + ((stopb-startb)/2)
        dist = centera-centerb
        proximities.append(dist)

    return proximities


do_overlap(options.a,options.b)
do_overlap(options.b,options.a)
proximities_a = do_proximity(options.a,options.b)
proximities_b = do_proximity(options.b,options.a)

fig = figure()
logged = False
normed = True
bins = arange(-1000000, 1000000,10000)
hist(proximities_a,bins=bins,alpha=0.5,log=logged,normed=normed)
hist(proximities_b,bins=bins,alpha=0.5,log=logged,normed=normed)
xlabel('distance between features (bp)')
ylabel('normalized count')
show()

