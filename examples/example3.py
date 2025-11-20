import os
import subprocess
import tempfile
from Bio.SeqIO import parse
from Bio.Seq import reverse_complement
from Bio.Align import PairwiseAligner, Alignment
import numpy as np

aligner = PairwiseAligner(scoring="blastn")


def permutate_trace(reference: str, sanger_trace: str) -> str:
    """Permutate a trace with respect to the reference using MARS"""
    # As an input for MARS, we need the reference + all traces
    # We include traces in both directions, since MARS does not handle
    # reverse complements - see https://github.com/lorrainea/MARS/issues/17#issuecomment-2598314356
    len_diff = len(reference) - len(sanger_trace)
    padded_trace = sanger_trace

    # Simple way to discriminate between Sanger / full sequence sequencing
    # If sanger, we pad with Ns to the length of the reference
    # If full sequence, we don't pad
    if len_diff > 0 and (len(sanger_trace) / len(reference) < 0.8):
        padded_trace = sanger_trace + len_diff * "N"

    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = os.path.join(tmpdir, "input.fa")
        with open(input_path, "w") as f:
            f.write(f">ref\n{reference}\n")
            f.write(f">trace\n{padded_trace}\n")

        output_path = os.path.join(tmpdir, "output.fa")
        result = subprocess.run(['mars', '-a', 'DNA', '-m', '0', '-i', input_path, '-o', output_path, '-q', '5', '-l', '20', '-P', '1'], capture_output=True, text=True)  # fmt: skip

        if result.returncode != 0:
            raise RuntimeError(f"MARS failed:\n{result.stderr}")

        # read permutated trace (second sequence in the FASTA file)
        ref, trace = parse(output_path, "fasta")
        return str(trace.seq)


plasmid_sequence = str(next(parse("pFA6a-kanMX6.gbk", "genbank")).seq)

# Trace 1 spans the origin and is in the reverse direction (insertions and substitutions highlighted in caps)
trace_1 = "tcgccgcagccgaacgaccgagcgcagcgAgtcagtgagcgaggaagcggaagagcgcccaatacgTTcaaaccgcctctccccgcgcgttggccgattcattaatgcaggttaacctggcttatcgaaattaatacgactcactatagggagaccggcagatcc"

# Trace 2 does not span the origin and is in the forward direction (insertions and substitutions highlighted in caps)
trace_2 = "cactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatacggttatccacagaatcaggggataacgcaggaaagaacatgtgagcaaaaggccagcaaaa"

# Trace 3 is a full plasmid sequence (we just rotate, reverse complement, and insert a few bases)
trace_3 = plasmid_sequence[:1000] + "AACCGTGCTG" + plasmid_sequence[1000:] + "CC"

pwas = []
for trace in [trace_1, trace_2, trace_3]:
    fwd = permutate_trace(plasmid_sequence, trace)
    rvs = permutate_trace(plasmid_sequence, reverse_complement(trace))

    unique_fwd_values = set(fwd)
    # Pairwise-align and keep the best alignment
    fwd_alignment = next(aligner.align(plasmid_sequence, fwd))
    rvs_alignment = next(aligner.align(plasmid_sequence, rvs))

    best_alignment = (
        fwd_alignment if fwd_alignment.score > rvs_alignment.score else rvs_alignment
    )

    pwas.append(best_alignment)

msa = Alignment.from_alignments_with_same_reference(pwas)
# Replace padding Ns with gaps and drop regions where every sequence is a gap
# introduced solely by the padding step
msa_array = np.asarray(msa, dtype=str)
msa_array[msa_array == "N"] = "-"
cols_to_keep = ~(np.all(msa_array == "-", axis=0))
if np.any(cols_to_keep):
    msa_array = msa_array[:, cols_to_keep]
msa = Alignment(["".join(row) for row in msa_array])

with open("example3-MSA.maf", "w") as f:
    f.write(msa.format("maf"))
