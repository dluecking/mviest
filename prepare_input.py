import os
import argparse
import sys
import shutil
from shutil import copy2, copyfileobj
import gzip

PARENT_DIRECTORY = "results/"

parser = argparse.ArgumentParser(description="This script prepares input for a full run of mv_signals")
parser.add_argument("sample", help = "name of the sample")

parser.add_argument("-vc", "--virome_contigs")
parser.add_argument("-mc", "--metagenome_contigs")

parser.add_argument("-mcgz", '--metagenome_contigs_gz')

parser.add_argument("-vr1", "--virome_reads_1", help = "forward reads of the virome")
parser.add_argument("-vr2", "--virome_reads_2", help = "reverse reads of the virome")

parser.add_argument("-mr1", "--metagenome_reads_1", help = "forward reads of the metagenome")
parser.add_argument("-mr2", "--metagenome_reads_2", help = "reverse reads of the metagenome")

parser.add_argument("-fo", "--force_overwrite", action="store_true", help="overwrite files, even if existing")
args = parser.parse_args()


# force overwrite mode?
if args.force_overwrite:
    FORCE_OVERWRITE = True 
else:
    FORCE_OVERWRITE = False

# basic structure
SAMPLE = args.sample
OUTDIR = os.path.join(PARENT_DIRECTORY, SAMPLE)

print(f' preparing {OUTDIR}')

if os.path.exists(OUTDIR):
    if FORCE_OVERWRITE:
        print(f"An output folder for {SAMPLE} exists already. Overwriting...")
        shutil.rmtree(OUTDIR)
    else:
        print(f"An output folder for {SAMPLE} exists already. Consider using the --force_overwrite flag.")
        sys.exit()

os.mkdir(OUTDIR)

READS_DIR = os.path.join(OUTDIR, "reads")
CONTIGS_DIR = os.path.join(OUTDIR, "contigs")

os.mkdir(READS_DIR)
os.mkdir(CONTIGS_DIR)

# move contig files (if existing)

def copyContigs(infile, outfile):
    if infile and os.path.exists(infile):
        copy2(infile, outfile)

copyContigs(infile=args.virome_contigs, 
            outfile=os.path.join(CONTIGS_DIR, "virome.contigs.fasta"))
copyContigs(infile=args.metagenome_contigs, 
            outfile=os.path.join(CONTIGS_DIR, "metagenome.contigs.fasta"))
copyContigs(infile=args.metagenome_contigs_gz, 
            outfile=os.path.join(CONTIGS_DIR, "metagenome.contigs.filtered.fasta.gz"))

# move virome read files (if existing) 
def copyReads(r1, r2, r1_out, r2_out):
    # read 1
    if r1 and os.path.exists(r1):
        if ".gz" in r1:
            copy2(r1, r1_out)
        else:
            print(r1)
    # read 2       
    if r2 and os.path.exists(r2):
        if ".gz" in r2:
            copy2(r2, r2_out)
        else:
            print(r2)
    else:
        print(r2)
        
        
copyReads(r1=args.virome_reads_1, 
          r2=args.virome_reads_2,
          r1_out=os.path.join(READS_DIR, "virome.reads1.fastq.gz"),
          r2_out=os.path.join(READS_DIR, "virome.reads2.fastq.gz"))

copyReads(r1=args.metagenome_reads_1, 
          r2=args.metagenome_reads_2,
          r1_out=os.path.join(READS_DIR, "metagenome.reads1.fastq.gz"),
          r2_out=os.path.join(READS_DIR, "metagenome.reads2.fastq.gz"))

    

