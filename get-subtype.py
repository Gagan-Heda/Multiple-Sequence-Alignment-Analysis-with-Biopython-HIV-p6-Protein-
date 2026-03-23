#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

def hamming_distance(seq1, seq2):
    return sum(a != b for a, b in zip(seq1, seq2))

def main():
    if len(sys.argv) != 3:
        print("Usage: ./get-subtype.py <msa.fasta> <input_sequence>")
        sys.exit(1)

    infile, query = sys.argv[1], sys.argv[2]
    query_seq = Seq(query)
    subtype_dists = defaultdict(list)

    for rec in SeqIO.parse(infile, "fasta"):
        # Example: "B.FR.83.HXB2_LAI_IIIB_BRU_K03455" → split('.') → ["B", "FR", "83", "HXB2_LAI_IIIB_BRU_K03455"] → [0] → "B"
        # Calculates Hamming distance between query and this sequence
        # Appends the distance to the list for that subtype
        subtype = rec.id.split('.')[0]  # subtype = first token (e.g. 'B')  
        dist = hamming_distance(str(rec.seq), str(query_seq))
        subtype_dists[subtype].append(dist)

    # Calculates average distance for each subtype using a dictionary comprehension
    # generates dictionary as subtype and distance
    avg_dists = {st: sum(vals)/len(vals) for st, vals in subtype_dists.items()}

    #Finds the subtype with the smallest average distance
    best_subtype = min(avg_dists.keys(), key=lambda k: avg_dists[k])
    print(best_subtype)

if __name__ == "__main__":
    main()