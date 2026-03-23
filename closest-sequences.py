#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def hamming_distance(seq1, seq2):
    return sum(a != b for a, b in zip(seq1, seq2))

def main():
    if len(sys.argv) != 3:
        print("Usage: ./closest-sequences.py <msa.fasta> <input_sequence>")
        sys.exit(1)

    infile, query = sys.argv[1], sys.argv[2]
    query_seq = Seq(query)
    records = list(SeqIO.parse(infile, "fasta"))

    distances = []
    for rec in records:
        seq_str = str(rec.seq)
        dist = hamming_distance(seq_str, str(query_seq))
        distances.append((rec, dist))

    distances.sort(key=lambda x: x[1])

    for rec, dist in distances[:3]:
        print(rec.seq)

if __name__ == "__main__":
    main()