from Bio import SeqIO
import sys

if len(sys.argv) != 3:
    sys.exit(1)

input_fasta = sys.argv[1]
epitope = sys.argv[2]
#Parse FASTA file
records = list(SeqIO.parse(input_fasta, "fasta"))

match_count = 0
# Loop through sequences and check for epitope
for rec in records:
    seq = str(rec.seq)
    if epitope in seq:
        match_count += 1

fraction = match_count / len(records) if len(records) > 0 else 0

print(fraction)
