from Bio.Align import PairwiseAligner, Alignment
import numpy as np

# Basic usage with strings
reference_str = "ACGT"
seq1_str = "ACT"
seq2_str = "ACGGT"
seq3_str = "AT"

aligner = PairwiseAligner()
pwa = next(aligner.align(reference_str, seq1_str))

coords = np.array([
    [0, 1, 2, 3, 3, 4],
    [0, 1, 2, 3, 4, 5],
    [0, 1, 1, 1, 1, 2]
])
not_pwa = Alignment([reference_str, seq2_str, seq3_str], coords)

print("Pairwise alignment:")
print(pwa)

print("\nThree-sequence alignment:")
print(not_pwa)

# multiple sequence alignment
msa = Alignment.from_alignments_with_same_reference([pwa, not_pwa])
print("\nCombined multiple sequence alignment:")
print(msa)
