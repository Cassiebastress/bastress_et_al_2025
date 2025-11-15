from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner, Alignment

# Using SeqRecord objects with metadata
reference_seqr = SeqRecord(Seq("ACGT"), id="reference", description="desc 1")
seq1 = SeqRecord(Seq("ACGGT"), id="seq1", description="desc 2")
seq2 = SeqRecord(Seq("AT"), id="seq2", description="desc 3")

aligner = PairwiseAligner()
pwa1 = next(aligner.align(reference_seqr, seq1))
pwa2 = next(aligner.align(reference_seqr, seq2))

msa = Alignment.from_alignments_with_same_reference([pwa1, pwa2])

print("Multiple sequence alignment with metadata preserved\n", msa)
