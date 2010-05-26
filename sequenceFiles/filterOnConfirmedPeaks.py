#!/usr/bin/python

import optparse
import sys
import numpy as np
from bedparser import bedfile
"""
This script takes a small set of confirmed peaks (``--bed-confirmed``) and
estimates the average peak height of those peaks in ``--wig-confirmed``.
This is the threshold to be used in the second step below.

It then takes a set of candidate peaks (``--bed-candidate``) and the values
at each of these regions (``--wig-candidates``), and returns a bed file
(``--outbed``) of only those candidate peaks that corresponded to 
peaks equal to or greater than the threshold peak height determined in 
the first step.

Use the Genome Browser for intersecting bed files with wig files.
Note that the bed files are required to correctly parse the 
variableStep wig files.
"""

op = optparse.OptionParser()

op.add_option('--wig-confirmed',
              help="""A wig file that has been intersected with a 
              bed file of confirmed peaks"""
              ,dest='wigconfirmed')
op.add_option('--wig-candidate',
              help="""A wig file that has been intersected
              with a bed file of candidate peaks""",
              dest='wigcandidate')
op.add_option('--bed-confirmed',
              help="""The bed file of confirmed peaks which was 
              intersected with the full wig data set""",
              dest='bedconfirmed')
op.add_option('--bed-candidate',
              help="""The bed file of candidate peaks which was
              intersected with the full wig data set""",
              dest='bedcandidate')
op.add_option('--trackname',
              help='''track name for output bed file.  
              Be sure to escape any spaces.''',
              dest='trackname', 
              default='User track')
op.add_option('--outbed', help='output bed file',dest='outbed')
options,args = op.parse_args()

def extract_chr(x):
    'returns chromosome name from a wig file track line definition'
    return x.split()[1].split('=')[1]

def get_min_dist(fn):
    '''Returns the minimum distance between bed features'''
    bedconfirmed = bedfile(fn)
    pos = 0
    mindist = 1e35
    chr = None
    for b in bedconfirmed:
        if b.chr == chr:
            dist = b.start-pos
            if dist < mindist:
                mindist = dist
        else:
            chr = b.chr
        pos = b.stop
    return mindist

# command line parsing
fail = False
for option in op.option_list:
    od = option.dest
    if od is not None:
        if getattr(options,od, None) is None:
            print 'Must specify %s' % od
            fail = True
if fail:
    op.print_help()
    sys.exit()

# get minimum distance in bed file of confirmed peaks
mindist_confirmed = get_min_dist(options.bedconfirmed)

# setup
last_pos = 0
regionlist = []  # list of lists, each contains wig values from bed region
thisregion = []  # the latest region, reset upon reaching a new bed region

for line in open(options.wigconfirmed):
    L = line.rstrip().split('\t')
    try:
        pos = int(L[0])
        value = float(L[1])
    except ValueError: # wasn't a normal data line
        continue
    if pos - last_pos >= mindist_confirmed: # reached a new region
        if len(thisregion) > 0: # first time through will be an empty list
            regionlist.append(thisregion)
            thisregion = [] # reset
    thisregion.append(value) # always append the value
    last_pos = pos

# done . . . but need to add the last region to the list.
regionlist.append(thisregion)

maxlist = [max(i) for i in regionlist]
threshold = np.median(maxlist)
print 'mean peak height: %s' % np.mean(maxlist)
print 'median peak height: %s' % np.median(maxlist)

# minimum distance between features in bed file of candidate regions
mindist_candidate = get_min_dist(options.bedcandidate)

# setup
last_pos = 0
regionlist = []
thisregion = []
thesevalues = [] 
chrs = [] # bed file output needs chr on each line
last_chr = None

for line in open(options.wigcandidate):
    L = line.rstrip().split('\t')
    try:
        pos = int(L[0])
        value = float(L[1])
    except ValueError:
        if 'chrom=' in line:
            chr = extract_chr(line)
        continue 

    # check to see if this position is very far from the last one.
    if (pos - last_pos <= mindist_candidate) and (chr == last_chr):
        thisregion.append(pos)
        thesevalues.append(value)
    else: # you are now in a new region 
        if len(thisregion) > 0:
            if max(thesevalues) >= threshold or chr != last_chr:
                regionlist.append(thisregion)
                chrs.append(chr)
        thisregion = []  # reset if not close enough to last one.
        thesevalues = [] # built in parallel with thisregion
    last_pos = pos
    last_chr = chr

# export to bed file.
fout = open(options.outbed,'w')
fout.write('track name="%s"\n' % options.trackname)
for chr,region in zip(chrs,regionlist):
    fout.write('%s\t%s\t%s\n' %(chr,region[0],region[-1]))
fout.close()
