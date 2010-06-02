#!/usr/bin/python
"""
Script to clean a GFF file downloaded from FlyBase.  Removes the sequence
and any features that are one of the types listed below.
"""

import sys
import os
import optparse
import GFFutils
usage = """
Cleans a GFF file so that it does not have "orthologous_to", "pcr_product", or
"BAC_cloned_genomic_insert" features.

Optionally adds a "chr" to the beginning of each chromosome.

Writes to stdout.
"""
op = optparse.OptionParser(usage=usage)
op.add_option('--addchr', action='store_true', help='Prefix each chromosome with "chr" in the output file.')

options,args = op.parse_args()

gfffn = args[0]

# this is a list of featuretypes that you want to remove.
culled_features = ['orthologous_to', 'pcr_product', 'BAC_cloned_genomic_insert']

f = GFFutils.GFFFile(gfffn)
for feature in f:
    if feature.start is None:
        continue
    if feature.stop is None:
        continue
    if feature.start > feature.stop:
        continue
    if feature.featuretype in culled_features:
        continue
    if options.addchr:
        feature.chr = 'chr'+feature.chr
    sys.stdout.write(feature.tostring())
   
