from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner, Alignment

# Reference and query sequences with metadata
reference = SeqRecord(Seq("ATGCCTA"), id="reference", description="desc1")
seq1 = SeqRecord(Seq("ATGCTAGCTA"), id="seq1", description="desc2")
seq2 = SeqRecord(Seq("ATTA"), id="seq2", description="desc3")

# Generate PWAs
aligner = PairwiseAligner()
pwa1 = next(aligner.align(reference, seq1))
pwa2 = next(aligner.align(reference, seq2))

# Display PWAs:
print("PWA 1:")
print(pwa1[0])
print(pwa1[1])

print("\nPWA 2:")
print(pwa2[0])
print(pwa2[1])

# Create and display the MSA
msa = Alignment.from_alignments_with_same_reference([pwa1, pwa2])
print("\nMSA:")
print(msa.format("fasta"))
