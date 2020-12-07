from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

# for each SV
orecs = []
for ff in snakemake.input:
    # read fasta file and keep largest contig
    contig = ''
    for record in SeqIO.parse(ff, "fasta"):
        if len(record.seq) > len(contig):
            contig = str(record.seq)
    # save for output fasta
    if contig != '':
        orecs.append(SeqRecord(Seq(contig),
                               id=ff, description=''))

# write fasta file with those reads
SeqIO.write(orecs, snakemake.output[0], "fasta")
