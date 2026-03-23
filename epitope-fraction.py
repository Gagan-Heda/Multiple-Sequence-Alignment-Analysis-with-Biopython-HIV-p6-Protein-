#!/usr/bin/env python3
import sys
from Bio import SeqIO

def main():
    if len(sys.argv) != 3:
        print("Usage: ./epitope-fraction.py <msa.fasta> <epitope>")
        sys.exit(1)

    infile, epitope = sys.argv[1], sys.argv[2]
    records = list(SeqIO.parse(infile, "fasta"))

    total = len(records)
    #Counts how many sequences contain the epitope
    # For each record, converts its sequence to a string:
    # Check if the given epitope is present or not
    # Then take sum for present epitope(1 if true)
    matches = sum(1 for rec in records if epitope in str(rec.seq))
    print(matches / total)

if __name__ == "__main__":
    main()