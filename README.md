# Multiple Sequence Alignment Analysis with Biopython (HIV p6 Protein)

Python-based workflow for analyzing multiple sequence alignments (MSA) of HIV p6 protein sequences. Includes gap filtering, consensus generation, epitope analysis, and sequence similarity-based classification.

---

## Overview

HIV p6 protein exhibits high sequence variability and contributes to viral pathogenicity. This project applies computational analysis to an MSA of p6 sequences to extract biologically meaningful insights.

**Goals:**

- Filter highly gapped alignment regions  
- Generate consensus sequences  
- Quantify epitope presence  
- Identify sequences most similar to a query  
- Classify sequences by subtype similarity  

---

## Scripts

- `compress-gaps.py`  
  - Removes columns where ≥67% of sequences are gaps  
  - Produces filtered MSA for downstream analysis  

- `print-consensus.py`  
  - Computes consensus sequence across alignment  
  - Selects most frequent amino acid per position  

- `epitope-fraction.py`  
  - Calculates fraction of sequences containing a given epitope  
  - Treats gaps as mismatches  

- `closest-sequences.py`  
  - Identifies top 3 most similar sequences to a query  
  - Uses Hamming distance as similarity metric  

- `get-subtype.py`  
  - Predicts subtype of a query sequence  
  - Computes average distance to sequences grouped by subtype  
  - Assigns subtype with smallest average distance  

---

## Methodology

### Alignment Processing
- Read FASTA-formatted MSA  
- Filter columns based on gap frequency  
- Preserve sequence integrity after filtering  

### Consensus Calculation
- Evaluate each alignment position  
- Select most frequent amino acid (including gaps if present)  

### Epitope Analysis
- Search for exact subsequence matches  
- Compute fraction of sequences containing the epitope  

### Sequence Similarity
- Compute Hamming distance between sequences  
- Rank sequences by similarity to query  

### Subtype Classification
- Extract subtype from sequence identifiers  
- Compute average distance between query and subtype groups  
- Assign subtype with minimum average distance  

---

## Technologies

- Python  
- Biopython  
- NumPy  

---

## Input

- MSA in FASTA format  
- Query sequences and epitopes via command-line arguments  

---

## Output

- Filtered MSA files  
- Consensus sequence  
- Epitope frequency values  
- Ranked similar sequences  
- Predicted subtype classifications  

---

## Key Features

- Automated MSA preprocessing  
- Consensus sequence computation  
- Epitope-level frequency analysis  
- Sequence similarity ranking  
- Subtype classification based on distance metrics  

---

## Skills Demonstrated

- Sequence analysis and bioinformatics workflows  
- Multiple sequence alignment processing  
- String and sequence manipulation  
- Distance-based similarity metrics  
- Classification using biological metadata  
- Python scripting with Biopython  

---

## Notes

- Gaps are treated as mismatches in comparisons  
- Hamming distance assumes aligned sequences of equal length  
- Subtype classification uses a simple distance-based heuristic  

---

## Author

**Gagan Heda**
