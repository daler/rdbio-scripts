#!/usr/bin/python

import optparse
import sys

__doc__ ="""
Script to convert a GFF file (e.g., from Nimblegen) into a bedGraph format that
can be used in the Genome Browser. """

op = optparse.OptionParser()
op.add_option('-i',dest='infn',help='Input GFF file')
op.add_option('--trackname', dest='trackname', help='Optional track name')
op.add_option('-o',dest='outfn',help='Output bedgraph file')
__doc__ += op.format_help()

def gff2bedgraph(infn, trackname=None, outfn=None):
    """Converts a GFF file from Nimblegen into a bedGraph, 
    wig-like file that can be viewed in the genome browser.

    *infn*
            
        input GFF file

    *outfn*
        
        Output bedGraph file

    *trackname*

        Optional track name
    """ 
    if outfn is not None:
        fout = open(outfn,'w')
    else:
        fout = sys.stdout
    

    if trackname is not None:
        fout.write('track type=bedGraph ')
        fout.write('name=%s'% trackname)
        fout.write('\n')
        
    for line in open(infn):
        L = line.rstrip().split('\t')
        chr,prog,id, start,stop,value,_,_,probe = L
        fout.write('\t'.join(  [chr, start,stop,value]) + '\n')
    
    if outfn is not None:
        fout.close()

if __name__ == "__main__":
    
    options,args = op.parse_args()
    gff2bedgraph(options.infn, options.trackname)
