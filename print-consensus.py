#!/usr/bin/env python3
import sys
from Bio import AlignIO
from collections import Counter
from Bio.Seq import Seq

def main():
    if len(sys.argv) != 2:
        print("Usage: ./print-consensus.py <msa.fasta>")
        sys.exit(1)

    infile = sys.argv[1]
    alignment = AlignIO.read(infile, "fasta")
    aln_len = alignment.get_alignment_length()

    consensus_chars = []
    for i in range(aln_len): #looping through each position in the alignment 
        col = [str(rec.seq[i]) for rec in alignment] #extracts the character at position i
        #Finds the most common character at this position 
        # Creates a Counter object that counts frequencies of each character 
        # Returns a list of ingle most common characters,count pair
        consensus_chars.append(Counter(col).most_common(1)[0][0])

    #Creates a BioPython Seq object by joining all consensus characters into a single string.
    consensus_seq = Seq(''.join(consensus_chars))
    print(consensus_seq)

if __name__ == "__main__":
    main()