from Bio import SeqIO
from collections import Counter
import sys
import numpy as np
#It requires 1 argument
if len(sys.argv) != 2:
    sys.exit(1)

input_fasta = sys.argv[1]

#Parse the FASTA alignment
records = list(SeqIO.parse(input_fasta, "fasta"))
seqs = [str(r.seq) for r in records]
#Convert alignment into a NumPy array
msa_array = np.array([list(s) for s in seqs])

#Loop through each alignment column
consensus = ""
for i in range(msa_array.shape[1]):
    #column is extracted
    col = msa_array[:, i]
    #Count characters
    counts = Counter(col)
    #Take the most common one
    consensus += counts.most_common(1)[0][0]
print(consensus)