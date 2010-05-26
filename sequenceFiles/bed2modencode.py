'''This script converts a BED-format file into one for modENCODE.

INPUT
-----
A BED file with tab-delimited syntax like this::

    browser position chrX:21,464,982-21,648,868
    browser hide all
    track name='piRNA (Yin and Lin)' itemRgb='On' visibility=4
    chr3R    5926673  5926694    piRNA(YL)    70    +    5926673    5926694    0,0,0
    chr3R    2147228  121472303  piRNA(YL)    66    +    21472281   21472303   0,0,0
    chr2R    9216924  9216946    piRNA(YL)    59    -    9216924    9216946    0,0,0
    chr3R    6233869  6233891    piRNA(YL)    39    +    6233869    6233891    0,0,0
    chr3L    10358380 10358402   piRNA(YL)    37    +    10358380   10358402   0,0,0
    chr3R    5916902  5916923    piRNA(YL)    32    +    5916902    5916923    0,0,0
    chr3R    16561698 16561719   piRNA(YL)    29    +    16561698   16561719   0,0,0

OUTPUT
------
A modENCODE-specific format that looks like this::

    glyph   = segments
    bgcolor = green
    key     = Mapped Expressed Tags
    
    mRNA hox1 Chr1:1..100 Type=UTR
    mRNA hox1 Chr1:101..200,300..400,500..800 Type=CDS
    mRNA hox1 Chr1:801..1000 Type=UTR

Columns are 
Feature type: any description
Feature name: any name, displayed when there's room for it
Position: chrX:1..100 or chrX:1-100.  I think you might have to remove 
          the "chr" part for modENCODE.
'''

import gzip
import sys

if 0:
    f = gzip.open('upload_to_UCSC.bed.gz')
    fout = open('upload_to_modENCODE.txt','w')
    for line in f:
        if 'browser' in line or 'track' in line:
            continue
        
        if 'small_RNA' not in line:
            continue

        if 'chrX' not in line:
            continue
        # chr3R    16561698 16561719   piRNA(YL)    29    +    16561698   16561719   0,0,0
       
        chrom,start,stop,desc,value,strand,start2,stop2,color = line.split('\t')
        chrom = chrom.replace('chr','')
        modENCODE_position = chrom + ':' + start + '-' + stop
        newline = '\t'.join(['small RNA',desc,modENCODE_position]) + '\n'
        fout.write(newline)      
    fout.close()
    f.close()

if 1:
    
    f = open('GSM322338.gff')
    fout = open('bed-from-encode.bed','w')

    for line in f:
        if line[0]=='#':
            continue
        if line[0:9] == 'Sequence:':  # not sure what these lines are for.
            continue
        line = line.rstrip()
        if len(line)< 1:
            continue
        chrom,desc,match,start,stop,value,strand,dot,data = line.split('\t')
        chrom = 'chr'+chrom
        chrom = chrom.replace('XHet','X')
        desc = desc.replace(' ','_')
        newline = '\t'.join( [chrom, start,stop,desc,value,strand,start,stop,'0,0,0'])+ '\n'

        fout.write(newline)

    f.close()
    fout.close()



