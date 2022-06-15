from Bio import SeqIO
import sys
import os
from numpy import append
import pandas as pd

# usage:
# python contig_selector.py {sample}

SAMPLE = sys.argv[1]
CONTIG_DIR = os.path.join("results/", SAMPLE, "contigs/")
ALL_CONTIGS_FILE = os.path.join(CONTIG_DIR, "virome.contigs.filtered.fasta")
SUMMARY_FILE = os.path.join("results/", SAMPLE, "contig_summary_"+SAMPLE+".tsv")
MV_CONTIG_FILE = os.path.join(CONTIG_DIR, "true.mvome.fasta")
VIROME_CONTIG_FILE = os.path.join(CONTIG_DIR, "true.virome.fasta")


input_handle = open(ALL_CONTIGS_FILE, "r")
records = SeqIO.parse(input_handle, "fasta")


df = pd.read_csv(SUMMARY_FILE)

df_mvome = df[df.final_label != "viral"]
df_virome = df[df.final_label == "viral"]

mv_records = []
viral_records = []


for record in records:
    if record.id in list(df_mvome.contig_id):
        mv_records.append(record)
    elif record.id in list(df_virome.contig_id):
        viral_records.append(record)
    else:
        print(f'whats up with {record.id}')
        pass

SeqIO.write(mv_records, MV_CONTIG_FILE, "fasta")
SeqIO.write(viral_records, VIROME_CONTIG_FILE, "fasta")
