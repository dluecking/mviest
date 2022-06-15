from Bio import SeqIO
import sys
import subprocess as sp
import os
import csv

contig_file = sys.argv[1]
out_file = sys.argv[2]
tmp_dir = sys.argv[3]

def is_circular(record):
    """Test if this contig can be circularized."""
    result = False
    contig_split_position = int(len(record.seq) / 2)

    contig_fragment_a_path = os.path.join(tmp_dir, "%s-a.fasta" % record.id)
    contig_fragment_a_seq = record.seq[:contig_split_position]
    with open(contig_fragment_a_path, mode='w') as fh:
        fh.write('>a\n')
        fh.write(str(contig_fragment_a_seq) + '\n ')

    contig_fragment_b_path = os.path.join(tmp_dir, "%s-b.fasta" % record.id)
    contig_fragment_b_seq = record.seq[contig_split_position:]
    with open(contig_fragment_b_path, mode='w') as fh:
        fh.write('>b\n')
        fh.write(str(contig_fragment_b_seq) + '\n ')

    cmd = [
        'nucmer',
        '-f',  # only forward strand
        '-l', '40',  # increase min match length to 40 bp
        '-p', record.id,
        "%s-a.fasta" % record.id,
        "%s-b.fasta" % record.id
    ]
    proc = sp.run(
        cmd, 
        cwd=str(tmp_dir), 
        stdout=sp.DEVNULL, 
        stderr=sp.DEVNULL, 
        universal_newlines=True
        )
    
    if(proc.returncode != 0):
        print("error")

    has_match = False
    with open(os.path.join(tmp_dir, "%s.delta" % record.id), "r") as fh:
        for line in fh:
            line = line.rstrip()
            print(line)
            if(line[0] == '>'):
                has_match = True
            elif(has_match):
                cols = line.split(' ')
                if(len(cols) == 7):
                    start_a = int(cols[0])
                    end_a = int(cols[1])
                    start_b = int(cols[2])
                    end_b = int(cols[3])
                    mismatches = int(cols[4])
                    alignment_a = end_a - start_a + 1
                    alignment_b = end_b - start_b + 1
                    if(alignment_a == alignment_b
                            and alignment_a > 50
                            and (mismatches / alignment_a) < 0.05
                            and end_b == len(contig_fragment_b_seq)
                            and start_a == 1):
                        result = True
                        break
    return(result)


handle = open(contig_file, "r")
records = SeqIO.parse(handle, "fasta")

# create tmp_dir
os.mkdir(tmp_dir)

id = []
circularity = []

for record in records:
    circularity.append(is_circular(record))
    id.append(record.id)
    for file in os.scandir(tmp_dir):
        if file.name.endswith(".fasta"):
            os.remove(file)

with open(out_file, "w") as f:
    writer = csv.writer(f)
    writer.writerows(zip(id, circularity))
