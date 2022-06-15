#Usage: python filter_contigs.py filename min_contig_length
from Bio import SeqIO
import sys
import os


CONTIG_FILE=sys.argv[1]
MIN_LEN=int(sys.argv[2])

handle = open(CONTIG_FILE, "r")
l = SeqIO.parse(handle, "fasta")

for s in l:
    if len(s.seq) >= MIN_LEN:
        print(">" + s.id)
        print(s.seq)
