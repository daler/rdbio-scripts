"""
All piRNAs except the 10% of reads corresponding to miRNAs, rRNAs, tRNAs, other
ncRNAs, and the sense strand of annotated genes were mapped to Release 5 and
the telomeric X-TAS repeat L03284.

Nucleotides corresponding to the 50 end of each piRNA were weighted according
to N/M with N = cloning frequency andM= number of genomic mappings. 

We used a 5 kb sliding window to identify all regions with densities greater
than 1 piRNA/kb. Windows within 20 kb of each other were collapsed into
clusters. Clusters with at least 5 piRNAs that uniquely matched to the cluster
were retained.

Supplemental:

For our analysis we exclusively used piRNAs matching the Release5 genome
assembly 100%. Excluded from our analysis was the “Uextra” file from Release5, which
includes short, un-assembled shotgun reads with low sequence quality and often
unverified origin. Less than 10% of the piRNAs that matched the Release5 genome
assembly uniquely had additional mappings in Uextra file, supporting the claim that
these sequences can be used to unambiguously identify the genomic origins of piRNAs.


So.

This implies:

1. reporting all multimappers
2. counting duplicates and storing their count in the BED score
3. counting different hits in the genome and storing that number . . . where? BED name?
4. intersecting with a GFF file and removing all genic transcripts.  With the FlyBase GFF, just use
   "gene" features.
5. tweaking the clustering code to handle BED scores

Hopefully 1-3 can be handled by SAMtools, and 4 with BEDtools.
"""


