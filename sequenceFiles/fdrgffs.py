#!/usr/bin/python

__doc__="""
Script to parse GFF files as output from NimbleGen's peak finder and write as a bed file.

If the feature has an FDR > ``--fdr``, then output that feature as a BED-format feature.

"""
import optparse, sys
op = optparse.OptionParser(usage='')
op.add_option('--fdr',dest='fdr',type=float,help='FDR, as a fraction')
op.add_option('-o', dest='outfn',help='Bed file to write out to')
op.add_option('-i', dest='infn',help='GFF file to parse')
__doc__ += op.format_help()

def fdrgffs(infn, outfn, fdr):
    """
    *infn*
        
        Input GFF file, from nimblegen

    *outfn*

        Output BED file

    *fdr*

        Float (0,1). Bed file will contain all features in the GFF with an FDR
        higher than *fdr*
    """
    def parsefdr(x):
        """
        From a GFF line, parse out the FDR.
        """
        L = line.rstrip().split('\t')
        details = L[-1].split(';')
        for d in details:
            if 'fdr' in d:
                fdr = d.split('=')[1]
                fdr = fdr.replace('attr_0','')
                return float(fdr)

    out = open(outfn,'w')
    for line in open(infn):
        if line.startswith('#'):
            continue
        if parsefdr(line) > fdr:
            continue
        L = line.strip().split('\t')
        chr = L[0]
        start = L[3]
        stop = L[4]
        out.write('%s\t%s\t%s\n' % (chr,start,stop))
    out.close()

if __name__ == "__main__":
    
    options,args = op.parse_args()
    fdrgffs(options.infn, options.outfn, options.fdr)

    
