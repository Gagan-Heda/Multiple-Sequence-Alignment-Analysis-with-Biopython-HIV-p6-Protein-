from Bio import SeqIO
import sys
import numpy as np

#It requires 2 arguments, 1 input + 1 output
if len(sys.argv) != 3:
    sys.exit(1)

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

#Parse the FASTA alignment
records = list(SeqIO.parse(input_fasta, "fasta"))

seqs = [str(r.seq) for r in records]
nseq = len(seqs)
seqlen = len(seqs[0])

#Convert the alignment into a NumPy matrix
msa_array = np.array([list(s) for s in seqs])
#Compute gap fraction for each column
gap_fraction = np.mean(msa_array == "-", axis=0)
#Decide which columns to keep
keep_cols = gap_fraction < 0.67
#Build compressed sequences
compressed_seqs = ["".join(msa_array[i, keep_cols]) for i in range(nseq)]

#Replace each SeqRecord’s sequence to compressed sequence
for i, rec in enumerate(records):
    rec.seq = rec.seq.__class__(compressed_seqs[i])

SeqIO.write(records, output_fasta, "fasta")

print(f"Compressed alignment saved to {output_fasta}")
print(f"Original length: {seqlen}; Compressed length: {sum(keep_cols)}")