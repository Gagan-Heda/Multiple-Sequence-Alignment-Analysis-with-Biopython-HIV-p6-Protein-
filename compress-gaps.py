#!/usr/bin/env python3
import sys
import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    if len(sys.argv) != 3:
        print("Usage: ./compress-gaps.py <input_msa.fasta> <output_msa.fasta>")
        sys.exit(1)

    infile, outfile = sys.argv[1], sys.argv[2]

    #Reads the input FASTA file as a multiple sequence alignment using BioPython's AlignIO.read() method.
    alignment = AlignIO.read(infile, "fasta")

    #Gets the number of sequences in the alignment
    nseq = len(alignment)
    #Gets the length of the alignment
    aln_len = alignment.get_alignment_length()

    # For each sequence record in the alignment, converts its sequence to a string and 
    # then to a list of characters 
    # Converts it into a numpy 2D array where rows = sequences and columns = alignment positions
    mat = np.array([list(str(rec.seq)) for rec in alignment])
    
    #This is a boolean array where the condition is true if cell is "-"
    # Then converts into numeric value and then computes mean
    gap_fraction = np.mean(mat == '-', axis=0)

    #Boolean array where columns whose gap fraction is less than 67% is kept
    keep_cols = gap_fraction < 0.67
    filtered_mat = mat[:, keep_cols]

    # Creates a new SeqRecord with the filtered sequence
    # keeping the original ID and description
    # Appends each new record to a list
    new_records = []
    for i, rec in enumerate(alignment):
        new_seq = Seq(''.join(filtered_mat[i]))
        new_records.append(SeqRecord(new_seq, id=rec.id, description=rec.description))

    SeqIO.write(new_records, outfile, "fasta")

if __name__ == "__main__":
    main()