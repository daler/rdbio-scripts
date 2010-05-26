#!/usr/bin/env python
import bedparser
import optparse

op = optparse.OptionParser()
op.add_option('-i',dest='input',help='input gff file')
op.add_option('-o',dest='output',help='output gff file')
options,args = op.parse_args()

fout = open(options.output,'w')
for feature in bedparser.gfffile(options.input):
    if feature.featuretype == 'gene':
        feature.attributes._attrs = ['ID']
    else:
        feature.attributes._attrs = ['ID','Parent']
    if not feature.chr.startswith('chr'):
        feature.chr = 'chr%s'%feature.chr
    fout.write(feature.tostring())

        
