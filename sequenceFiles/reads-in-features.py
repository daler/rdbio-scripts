#!/usr/bin/python

usage = r"""
Script that uses HTSeq to count reads in genes, introns, exons, etc.   

It will speed things up if you have a GFF file that only has gene, exon, and intron features.  For example,


    # Set up some environmental vars
    export GFF=your-gff.gff
    export OUT=smaller.gff

    # write gene features to new file $OUT
    awk -F "\t" '{if ($3 ~ /^gene/) print $0}' $GFF > $OUT

    # append exon features to $OUT
    awk -F "\t" '{if ($3 ~ /^exon/) print $0}' $GFF >> $OUT

    # append intron features to $OUT
    awk -F "\t" '{if ($3 ~ /^intron/) print $0}' $GFF >> $OUT


Then specify --gff=$OUT when calling this script.


Running this script results in a new file, `$OUTPREFIX.counts.report` which
contains the counts for the following regions: 

    * empty (i.e., reads not falling in any annotated regions)
    * genes
    * introns only
    * exons only
    * any exon
    * any intron
    * exons and introns (happens often with multiple isoforms)
    * spliced reads

In addition, if you use the --debug option, this script will also create WIG
files and spliced BED files for these regions.  Use with caution, because the
files can get quite large.

Expect the counting to take about a minute per 1M reads; allow extra time for
sorting and converting BED to WIG if --debug is specified.


Requirements:  bed2wig.py needs to be on your path
               HTSeq needs to be installed.
"""

import optparse
import sys
import HTSeq
import os
import time
import pdb

op = optparse.OptionParser(usage=usage)
op.add_option('--sam',help='Input SAM file (required)')
op.add_option('--gff',help='Input GFF file (required)')
op.add_option('--outprefix',help='Output BED file prefix for '
                              'bedfiles that will be created.  If a dir is in '
                              'the outprefix path, that dir will be created '
                              'if needed.  (required)')
op.add_option('--debug',action='store_true',
              help='Creates useful output, like BEDs and WIGs of counted reads, '
                   'useful for debugging or digging deeper into the returned counts. (optional)')
op.add_option('--verbose',action='store_true',help='Print progress to stderr (optional)')
op.add_option('--label',,help='Label for library that will be added to the top of count reports '
                              'and will be prefixed to track names if --debug is enabled (default '
                              'is to use the basename of the SAM file)')
options,args = op.parse_args()

reqs = ['sam','gff','outprefix']
for req in reqs:
    if getattr(options,req) is None:
        op.print_help()
        sys.stderr.write('\n\nRequired arg "%s" missing\n\n'%req)
        sys.exit(1)

if options.label is None:
    options.label = os.path.basename(options.sam)

outdir = os.path.split(options.outprefix)[0]
if not os.path.exists(outdir) and len(outdir) > 0:
    os.system('mkdir -p %s' % outdir)

class UnknownChrom( Exception ):
    pass
    

def sam2bed(r,splice=False):
    """
    Given a HTSeq.SAM_Alignment instance, returns a BED6 (if splice=False) or
    BED12 line (if splice=True).
    """
    chrom = r.iv.chrom
    start = r.iv.start
    stop = r.iv.end
    name = '.'
    score = 0
    strand = r.iv.strand
    if not splice:
        return '\t'.join(map(str,[chrom,start,stop,name,score,strand]))+'\n'
    else:
        itemRGB = '0,0,0'
        thickStart = r.iv.start
        thickEnd = r.iv.end
        blockCount = 0
        blockSizes = []
        blockStarts = []
        for co in r.cigar:
            if co.type == "M":
                blockCount += 1
                blockStarts.append(co.ref_iv.start-start)
                blockSizes.append(co.size)
        blockSizes = ','.join(map(str,blockSizes))
        blockStarts = ','.join(map(str,blockStarts))
        line = [chrom,start,stop,name,score,strand,thickStart,thickEnd,itemRGB,blockCount,blockSizes,blockStarts]
        line = map(str,line)
        return '\t'.join(line) + '\n'
    

# fail early on bad filenames
open(options.sam).close()

t0 = time.time()
counts = {}
output_beds = {}
featuretypes = ['gene','exon','intron','spliced','exon-and-intron','exon-only','intron-only','empty','total']

# Create new files for each class of features; write a header line too
for ft in featuretypes:
    counts[ft] = 0
    if options.debug:
        output_beds[ft] = open(options.outprefix+'.'+ft+'.bed','w')
        output_beds[ft].write('track name="%s reads"\n' % ft)


# Here we go: time to read in the GFF features as an HTSeq.GenomicArrayOfSets
featurecount = 0
gff = HTSeq.GFF_Reader(options.gff)
features = HTSeq.GenomicArrayOfSets([],stranded=False)
try:
    for f in gff:
        featurecount += 1
        if options.verbose:
            if featurecount % 5000:
                sys.stderr.write('\r%s GFF features imported...'%featurecount)
        if f.type not in featuretypes:
            continue

        # features.step_vectors is a dict of chroms . . . so if we don't
        # have an entry for this chrom yet, then add one.
        if f.iv.chrom not in features.step_vectors.keys():
            features.add_chrom(f.iv.chrom)

        # Add this feature's featuretype to the "features" GenomicArrayOfSets.
        # In other words, the interval spanned by this GFF feature is annotated
        # by the featuretype.
        try:
            features.add_value(f.type, f.iv)
        except IndexError:
            sys.stdout.write('%s has funky coords, skipping\n'%f)
            continue

# Deal with any GFF file reading errors
except ValueError as e:
    e.args += ( gff.get_line_number_string(), )
    raise

try:
    # Get the first read to see if we're dealing with paired-end data
    read_seq = HTSeq.SAM_Reader(options.sam)
    first_read = iter(read_seq).next()
    pe_mode = first_read.paired_end
    
    # Re-initialize read_seq depending on if it's paired-end data or not
    read_seq = HTSeq.SAM_Reader(options.sam)
    if pe_mode:
        read_seq = HTSeq.pair_SAM_alignments(read_seq)

    # Read counter, for feedback to user
    i = 0 
    total = 0
    # Here we go, through each read...
    for r in read_seq:
        spliced = False
        if not pe_mode:
            if not r.aligned:
                continue
            total += 1
            iv_seq = []

            # Check to see if it's spliced
            for co in r.cigar:
                if co.type == 'N':
                    spliced = True
                elif co.type == 'M':
                    iv_seq.append(co.ref_iv)
            
        else:
            if r[0] is not None and r[0].aligned:
                iv_seq = []
                for co in r[0].cigar:
                    if co.type == 'N':
                        spliced = True
                    elif co.type == 'M':
                        iv_seq.append(co.ref_iv)
            else:
                iv_seq = tuple()

            # TODO:  not sure if you need to count the mate as spliced as well
            # . . . so it's possible the splice counts are underestimating the
            # true number, for paired-end data
            if r[1] is not None and r[1].aligned:                
                iv_seq = itertools.chain( iv_seq, 
                    ( invert_strand( co.ref_iv ) for co in r[1].cigar if co.type == "M" ) )
            else:
                if ( r[0] is None ) or not ( r[0].aligned ):
                    continue            
        try:

            # This empty set will have GFF feature IDs "unioned" to it
            # below...
            fs = set()

            # Recall iv_seq is a tuple of intervals.  So iterate through 'em
            bed_intervals = []
            if len(iv_seq) == 0:
                print 'length 0 iv_seq'
            for iv in iv_seq:
                if iv.chrom not in features.step_vectors:
                    raise UnknownChrom

                # This gets the unique GFF features that are within the
                # interval of this CIGAR operation.  We're in a unique feature
                # (one that does not overlap any others) if the length of fs ==
                # 1.  For a full GFF with genes, introns, and exons, fs should
                # not be exactly 1, since if you're in a gene you're also
                # either in an exon or an intron
                for fs2 in features.get_steps( iv, values_only=True ):
                    fs = fs.union( fs2 )
                
                # Add the intervals to the ongoing list of bed features.
                bed_intervals.append('%s\t%s\t%s\t.\t0\t%s\n' % (iv.chrom,iv.start,iv.end,iv.strand))

            # We'll at least be counting the featuretypes found here that
            # overlap the read (e.g. intron, exon) ...but we're also interested
            # in some derived featuretypes.  Whether or not to add this read to
            # one of those derived featuretypes is determined below, where a
            # derived class is added to this list.
            featuretypes_to_count = list(fs)


            # DERIVED FEATURETYPES
            # spliced reads will also be counted in genes and exons.
            if spliced:
                featuretypes_to_count.append('spliced')
            
            # spliced reads skip the spliced part -- so they should not be included here.
            if ('exon' in fs) and ('intron' in fs): 
                featuretypes_to_count.append('exon-and-intron')

            # not sure why you'd have an exon outside of a gene, but oh well,
            # as long as it doesn't have 'intron' in it...
            if fs == set(['exon','gene']):
                featuretypes_to_count.append('exon-only')
            
            # same for introns
            if fs == set(['intron','gene']):
                featuretypes_to_count.append('intron-only')
            
            # nothing annotated here
            if 'gene' not in fs:
                featuretypes_to_count.append('empty')
            if 'empty' in featuretypes_to_count:
                if len(featuretypes_to_count) != 1:
                    if 'spliced' not in featuretypes_to_count:
                        pdb.set_trace()
            
            # increment all featuretypes that were found here
            for featuretype in featuretypes_to_count:
                counts[featuretype] += 1
                if featuretype == 'spliced':
                    splice=True
                else:
                    splice=False

                # Dispatch dictionary of files to write to, depending on
                # featuretype.  Converts to BED format on the fly, and tries to
                # be smart about splices.
                if options.debug:
                    output_beds[featuretype].writelines(sam2bed(r,splice))

            counts['total'] += 1

        except UnknownChrom:
            if not pe_mode:
                rr = r
            else: 
                rr = r[0] if r[0] is not None else r[1]
            if not quiet:
                sys.stderr.write( ( "Warning: Skipping read '%s', because chromosome " +
                    "'%s', to which it has been aligned, did not appear in the GFF file.\n" ) % 
                    ( rr.read.name, iv.chrom ) )

        i += 1
        if options.verbose:
            if i % 5000 == 0:
                sys.stderr.write('\r%d reads processed, %ds elapsed'%(i, time.time()-t0))
                sys.stderr.flush()

except ValueError as e:
    e.args += ( read_seq.get_line_number_string(), )
    raise
if options.verbose:
    sys.stderr.write('\n\n')

# write out counts to the report
fout = open(options.outprefix+'.counts.report','w')
fout.write(('%s'%options.label)+'\n')
for fn in sorted( counts.keys() ):
    fout.write("%s\t%d\n" % ( fn, counts[fn] ) )
fout.close()


if options.debug:
    if options.verbose:
        sys.stderr.write('Sorting output BED files and converting to WIG...\n')
    for featuretype,f in output_beds.items():
        fn = f.name
        if options.verbose:
            sys.stderr.write(fn+'\n')
            sys.stderr.flush()
        
        # make sure you flush the tempfile
        f.close()
        label = options.label
        os.system('sortBed -i %(fn)s > %(fn)s.sorted'%locals())
        os.system('bed2wig.py -i %(fn)s.sorted -o %(fn)s.wig --type bed --track="name=\"%(label)s-%(featuretype)s\""' % locals())
        os.unlink('%(fn)s.sorted'%locals())

