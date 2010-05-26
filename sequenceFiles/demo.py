fn = 'input.fa'

# method 1
seq = open(fn).read().strip() # remove any newlines at the end
k = 3
mer_dict = {}
for i in range(len(seq)-k):
    kmer = seq[i:i+k]
    count = mer_dict.setdefault(kmer,0)
    count += 1
    mer_dict[kmer] = count
kmers_found = mer_dict.keys()

# print output
kmers_found.sort()
for kmer in kmers_found:
    print '%s\t%s' % (kmer, mer_dict[kmer])


# method 2
#seq = open(fn).read().strip()
#k = 3
#mer_dict = {}
#fout = open('output.txt','w')
#for i in range(len(seq)):
#    kmer = seq[i:i+k]
#    fout.write(kmer+'\n')
#fout.close()
#os.system('sort | uniq -c output.txt > final.txt')

# Finding seqs in FASTQ file is similar . . .
# BioPython has a FASTQ parser:

from Bio import SeqIO
fastq_fn = 's_1_sequence.txt'
handle = open(fastq_fn)
parser = SeqIO.parse(handle,'fastq')
fastq_kmers = {}
for rec in parser:
    seq = rec.seq
    for i in range(len(seq)-k):
        kmer = rec.seq[i:i+k]
        count = fastq_kmers.setdefault(kmer,0)
        count += 1
        fastq_kmers[kmer] = count
handle.close()


# the fastq file might have some kmers that are not found in 
# the spike-in sequence; you can parse these out of the dictionaries:
kmers_in_fastq = fastq_kmers.keys()
unique_to_fastq = set(kmers_in_fastq).difference(set(kmers_found))
for i in unique_to_fastq:
    print '%s: %s' % (i,fastq_kmers[i])

