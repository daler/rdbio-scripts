#!/usr/bin/python

__doc__="""
Parse bed files containing multimapping and unique reads (in that
order!), removing any features from the multimappers that are found in the
uniques.  In doing so, adds colors to the BED files. Result is two new BED
files that contain mutually exclusive features.

"""

import bedparser
import sys, optparse, subprocess, tempfile, os,logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger('bedColorUniques')

op = optparse.OptionParser(usage='')
op.add_option('-a', dest='a', help='Input bed file 1: multireads')
op.add_option('-b', dest='b', help='Input bed file 2: unique reads')
op.add_option('--oa', dest='oa', help='Output bed file of colored multireads')
op.add_option('--ob', dest='ob', help='Output bed file of colored unique reads')
op.add_option('--uniquecolor', dest='uniquecolor', help='csv RGB tuple for uniques (e.g., 0,0,255)')
op.add_option('--multicolor', dest='multicolor', help='csv RGB tuple for multimappers (e.g., 128,0,0) ')
__doc__ += op.format_help()

def bedColorUniques(a,b,oa,ob,multicolor,uniquecolor):
    tmp1 = tempfile.mktemp()
    tmp2 = tempfile.mktemp()
    cmds = """awk -F "\\t"  'BEGIN {OFS="\\t"}{print $1, $2, $3, $4, $5, $6}' %s > %s""" %(a,tmp1)
    log.info(cmds)
    os.system(cmds)

    cmds = """awk -F "\\t"  'BEGIN {OFS="\\t"}{print $1, $2, $3, $4, $5, $6}' %s > %s""" %(b,tmp2)
    log.info(cmds)
    os.system(cmds)

    log.info('intersecting')
    tmp3 = tempfile.mktemp()
    cmds = """
    intersectBed \\
    -a %s \\
    -b %s \\
    -v \\
    -s \\
    > %s""" % (tmp1, tmp2, tmp3)
    os.system(cmds)
    
    log.info('using bedparser.py to recolor multireads') 
    outa = open(oa,'w')
    for i in bedparser.bedfile(tmp3):
        i.thickStart = i.start
        i.thickStop = i.stop
        i.itemRGB = multicolor
        outa.write(i.tostring())
    outa.close()
    
    log.info('using bedparser.py to recolor unique reads')
    outb = open(ob,'w')
    for i in bedparser.bedfile(b):
        i.thickStart = i.start
        i.thickStop = i.stop
        i.itemRGB = uniquecolor
        outb.write(i.tostring())
    outb.close()

if __name__ == "__main__":
    options,args = op.parse_args()
    bedColorUnique(options.a, options.b, options.oa, options.ob, options.multicolor, options.uniquecolor)

