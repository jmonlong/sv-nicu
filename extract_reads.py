import pysam
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

flank_size = 1000

# open bam file
bamfile = pysam.AlignmentFile(snakemake.input['bam'], "rb")

inb = open(snakemake.input['svs'], 'r')
outb = open(snakemake.output[0], 'w')

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
    print(line[0] + ':' + line[1] + '-' + line[2])
    # reads to extract
    reads_t = line[rnames_ii].split(',')
    print('\t' + str(len(reads_t)) + ' reads to extract')
    # fetch reads and look for target reads
    reads_reg = bamfile.fetch(reference=line[0], start=int(line[1]) - flank_size,
                              end=int(line[2]) + flank_size)
    recs = []
    for read in reads_reg:
        if read.query_name in reads_t:
            print('\t\t' + read.query_name)
            if read.query_alignment_sequence != None:
                print('\t\t' + str(len(read.query_alignment_sequence)) + 'bp')
                recs.append(SeqRecord(Seq(read.query_alignment_sequence),
                                      id=read.query_name, description=''))
    print('\t' + str(len(recs)) + ' reads extracted')
    # write fasta file with those reads
    out_fa = line[svid_ii] + '.reads.fa'
    SeqIO.write(recs, out_fa, "fasta")
    outb.write('\t'.join(line[:3]) + '\t' + out_fa + '\n')

inb.close()
outb.close()