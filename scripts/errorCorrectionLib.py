#!/usr/bin/env python
"""Library for error correction
"""
from serviceCourse.sequenceTools import reverse_complement
from serviceCourse.parsers import read_fasta


def make_degenerate_reference_iterator(input_sequence, block_size=1, step=6):
    """
    input_sequence: string, input nucleotide sequence
    out_path: string, path to directory to put new sequences with substituted degenerate characters
    block_size: not implemented
    step: number of bases between degenerate characters
    :return (subbed sequence, complement subbed sequence)
    """
    complement_sequence = reverse_complement(dna=input_sequence, reverse=False, complement=True)

    for s in xrange(0, step):
        positions = xrange(s, len(input_sequence), step)
        t_seq = list(input_sequence)
        c_seq = list(complement_sequence)
        for position in positions:
            t_seq[position] = "X"
            c_seq[position] = "X"
        yield ''.join(t_seq), ''.join(c_seq)


def write_degenerate_reference_set(input_fasta, out_path):
    # get the first sequence from the FASTA
    seq = ""
    for header, comment, sequence in read_fasta(input_fasta):
        seq += sequence
        break

    for i, s in enumerate(make_degenerate_reference_iterator(input_sequence=seq)):
        with open(out_path + "forward_sub{i}.txt".format(i=i), 'w') as f:
            f.write("{seq}".format(seq=s[0]))
        with open(out_path + "backward_sub{i}.txt".format(i=i), 'w') as f:
            f.write("{seq}".format(seq=s[1]))
