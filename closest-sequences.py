from Bio import SeqIO
import sys

#It requires 2 argument. 1 in fasta file and 2 is query sequence
if len(sys.argv) != 3:
    sys.exit(1)

input_fasta = sys.argv[1]
query_seq = sys.argv[2]
#Parse the FASTA alignment
records = list(SeqIO.parse(input_fasta, "fasta"))

#Computes the number of mismatched positions between two sequences
def hamming_distance(seq1, seq2):
    # Pad shorter sequence with gaps so lengths match
    length = min(len(seq1), len(seq2))
    return sum(a != b for a, b in zip(seq1[:length], seq2[:length])) + abs(len(seq1) - len(seq2))

distances = []
for rec in records:
    #sequence is converted to string
    seq = str(rec.seq)
    #Compute distance between query and each FASTA sequence
    dist = hamming_distance(seq, query_seq)
    distances.append((seq, dist))

# Sorts the sequence similar to lease similar
distances.sort(key=lambda x: x[1])

#Top 3 similar seq
for seq, dist in distances[:3]:
    print(seq)
