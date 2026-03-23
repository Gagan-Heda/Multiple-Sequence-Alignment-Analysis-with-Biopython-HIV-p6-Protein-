from Bio import SeqIO
import sys
from collections import defaultdict

#It requires 2 argument. 1 in fasta file and 2 is query sequence
if len(sys.argv) != 3:
    sys.exit(1)

input_fasta = sys.argv[1]
query_seq = sys.argv[2]

def hamming_distance(seq1, seq2):
    length = min(len(seq1), len(seq2))
    return sum(a != b for a, b in zip(seq1[:length], seq2[:length])) + abs(len(seq1) - len(seq2))

#Create dictionary to store sequences by subtype
subtype_sequences = defaultdict(list)

#Parse FASTA and group sequences by subtype
for rec in SeqIO.parse(input_fasta, "fasta"):
    header = rec.id 
    subtype = header.split('.')[0]
    subtype_sequences[subtype].append(str(rec.seq))

#Compute average distance for each subtype
avg_distances = {}
for subtype, seqs in subtype_sequences.items():
    distances = [hamming_distance(seq, query_seq) for seq in seqs]
    avg_distances[subtype] = sum(distances) / len(distances)

#Pick subtype with smallest average distance
best_subtype = min(avg_distances, key=avg_distances.get)
print(best_subtype)