#!/usr/bin/env python
"""Takes alignments from signalAlign and calls methylation status"""
from __future__ import print_function, division
import glob
import os
import sys
from argparse import ArgumentParser
from alignmentAnalysisLib import CallMethylation
from serviceCourse.parsers import read_fasta
from random import shuffle
from multiprocessing import Process, current_process, Manager

def parse_args():
    parser = ArgumentParser(description=__doc__)

    # query files
    parser.add_argument('--C_alignments', '-C', action='store',
                        dest='C_alns', required=False, type=str, default=None,
                        help="C files")
    parser.add_argument('--mC_alignments', '-mC', action='store',
                        dest='mC_alns', required=False, type=str, default=None,
                        help="mC files")
    parser.add_argument('--hmC_alignments', '-hmC', action='store',
                        dest='hmC_alns', required=False, type=str, default=None,
                        help="hmC files")
    parser.add_argument('--ref', '-r', required=True, action='store', type=str, dest='ref',
                        help="path to fasta reference file")
    parser.add_argument('-n', required=False, action='store', type=int, dest='n', default=100,
                        help='Max number of alignments from each category to look at')
    #parser.add_argument('--out', '-o', action='store', type=str, required=True, dest='out_file')

    return parser.parse_args()


def randomly_select_alignments(path_to_alignments, max):
    if path_to_alignments is None:
        return None
    else:
        alignments = [x for x in glob.glob(path_to_alignments) if os.stat(x).st_size != 0]
        shuffle(alignments)
        if len(alignments) < max:
            return alignments
        else:
            return alignments[:max]


def get_forward_mask(list_of_alignments):
    mask = []
    for alignment in list_of_alignments:
        if alignment.endswith(".backward.tsv"):
            mask.append(False)
        else:
            mask.append(True)
    return mask


def get_reference_sequence(path_to_fasta):
    seqs = []

    for header, comment, sequence in read_fasta(path_to_fasta):
        seqs.append(sequence)

    assert len(seqs) > 0, "Didn't find any sequences in the reference file"

    if len(seqs) > 1:
        print("[NOTICE] Found more than one sequence in the reference file, using the first one")

    return seqs[0]


def get_alignments_labels_and_mask(path_to_alignments, label, max):
    if path_to_alignments is None:
        return [], [], []
    else:
        alignments = randomly_select_alignments(path_to_alignments, max)
        labels = [label] * len(alignments)
        mask = get_forward_mask(alignments)
        return alignments, labels, mask


def main(args):
    args = parse_args()
    reference_sequence = get_reference_sequence(args.ref)

    C_alns, C_labels, C_mask = get_alignments_labels_and_mask(args.C_alns, "C", args.n)
    mC_alns, mC_labels, mC_mask = get_alignments_labels_and_mask(args.mC_alns, "E", args.n)
    hmC_alns, hmC_labels, hmC_mask = get_alignments_labels_and_mask(args.hmC_alns, "O", args.n)

    all_alignments = C_alns + mC_alns + hmC_alns
    all_labels = C_labels + mC_labels + hmC_labels
    all_masks = C_mask + mC_mask + hmC_mask

    correct_template_calls = 0
    correct_complement_calls = 0
    total_template_calls = 0
    total_complement_calls = 0

    for aln, label, forward_bool in zip(all_alignments, all_labels, all_masks):
        call_methyl_args = {
            "sequence": reference_sequence,
            "alignment_file": aln,
            "forward": forward_bool,
            "test": True,
            "label": label
        }
        c = CallMethylation(**call_methyl_args)
        correct_calls = c.test()
        correct_template_calls += correct_calls["template"]
        correct_complement_calls += correct_calls["complement"]
        total_template_calls += len(c.template_calls)
        total_complement_calls += len(c.complement_calls)

    template_accuracy = correct_template_calls / total_template_calls * 100
    complement_accuracy = correct_complement_calls / total_complement_calls * 100
    print("Template accuracy {t} Complement accuracy {c} {nt} template cytosines tested "
          "{nc} complement cytosines tested".format(t=template_accuracy, c=complement_accuracy,
                                                    nt=total_template_calls, nc=total_complement_calls),
          file=sys.stdout)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
