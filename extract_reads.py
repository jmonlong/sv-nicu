import argparse
import pysam
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

flank_size = 1000

parser = argparse.ArgumentParser(description='Extract reads for a chr')
parser.add_argument('-b', help='BAM file', required=True)
parser.add_argument('-v', help='TSV file with variants', required=True)
parser.add_argument('-c', help='chr name', required=True)
parser.add_argument('-o', help='output file', required=True)
args = parser.parse_args()

# open bam file
bamfile = pysam.AlignmentFile(args.b, "rb")

outb = open(args.o, 'w')
inb = open(args.v, 'r')

# headers
heads = next(inb)
heads = heads.rstrip().split()
rnames_ii = 0
svid_ii = 0
for ii in range(len(heads)):
    if heads[ii] == 'RNAMES':
        rnames_ii = ii
    if heads[ii] == 'svid':
        svid_ii = ii
outb.write('chr\tstart\tend\treads\n')

# for each SV
for line in inb:
    line = line.rstrip().split('\t')
    # skip if not the specified chromosome
    if line[0] != args.c:
        continue
    print(line[0] + ':' + line[1] + '-' + line[2])
    # reads to extract
    reads_t = line[rnames_ii].split(',')
    print('\t' + str(len(reads_t)) + ' reads to extract')
    # fetch reads and look for target reads
    reads_reg = bamfile.fetch(reference=line[0],
                              start=int(line[1]) - flank_size,
                              end=int(line[2]) + flank_size)
    recs = []
    for read in reads_reg:
        if read.query_name in reads_t:
            print('\t\t' + read.query_name)
            if read.query_alignment_sequence is not None:
                print('\t\t' + str(len(read.query_alignment_sequence)) + 'bp')
                recs.append(SeqRecord(Seq(read.query_alignment_sequence),
                                      id=read.query_name, description=''))
    print('\t' + str(len(recs)) + ' reads extracted')
    # write fasta file with those reads
    if len(recs) > 0:
        out_fa = args.c + '_' + line[svid_ii] + '.reads.fa'
        SeqIO.write(recs, out_fa, "fasta")
        outb.write('\t'.join(line[:3]) + '\t' + out_fa + '\n')

inb.close()
outb.close()
