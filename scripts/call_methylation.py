#!/usr/bin/env python
"""Takes alignments from signalAlign and calls methylation status"""
from __future__ import print_function, division
import glob
import os
import sys
import numpy as np
from argparse import ArgumentParser
from alignmentAnalysisLib import CallMethylation
from serviceCourse.parsers import read_fasta
from random import shuffle
from multiprocessing import Process, current_process, Manager


def parse_args():
    parser = ArgumentParser(description=__doc__)

    # query files
    parser.add_argument('--input', '-i', action='store',
                        dest='in_files', required=False, type=str, default=None,
                        help="files, with file-suffix")
    parser.add_argument('--ref', '-r', required=False, action='store', type=str, dest='ref', default=None,
                        help="path to fasta reference file")
    parser.add_argument('--positions', '-p', required=False, action='store', type=str, dest='positions',
                        help='positions file')
    parser.add_argument('-n', required=False, action='store', type=int, dest='n', default=100,
                        help='Max number of alignments from each category to look at')
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=4, type=int, help="number of jobs to run concurrently")
    parser.add_argument('--out', '-o', action='store', type=str, required=True, dest='out')

    return parser.parse_args()


def randomly_select_alignments(path_to_alignments, max):
    alignments = [x for x in glob.glob(path_to_alignments) if os.stat(x).st_size != 0]
    if len(alignments) == 0:
        print("[error] Didn't find any alignment files here {}".format(path_to_alignments))
        sys.exit(1)

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


def get_alignments_labels_and_mask(path_to_alignments, max):
    alignments = randomly_select_alignments(path_to_alignments, max)
    mask = get_forward_mask(alignments)
    return alignments, mask


def parse_substitution_file(substitution_file):
    fH = open(substitution_file, 'r')
    line = fH.readline().split()
    forward_sub = line[0]
    forward_pos = map(np.int64, line[1:])
    line = fH.readline().split()
    backward_sub = line[0]
    backward_pos = map(np.int64, line[1:])
    return (forward_sub, forward_pos), (backward_sub, backward_pos)


def run_methyl_caller(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            c = CallMethylation(**f)
            c.write()
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def main(args):
    args = parse_args()

    if args.ref is not None:
        reference_sequence = get_reference_sequence(args.ref)
    else:
        reference_sequence = None

    alns, forward_mask = get_alignments_labels_and_mask(args.in_files, args.n)

    out_file = args.out

    if args.positions is not None:
        positions = {}
        f, b = parse_substitution_file(args.positions)
        positions['forward'] = f[1]
        positions['backward'] = b[1]
    else:
        positions = None

    workers = args.nb_jobs
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()
    jobs = []

    for aln, forward_bool in zip(alns, forward_mask):
        call_methyl_args = {
            "sequence": reference_sequence,
            "alignment_file": aln,
            "forward": forward_bool,
            "out_file": out_file,
            "positions": positions,
        }
        work_queue.put(call_methyl_args)
        #c = CallMethylation(**call_methyl_args)
        #correct_calls = c.test()
        #c.write(out_file=out_file)
        #correct_template_calls += correct_calls["template"]
        #correct_complement_calls += correct_calls["complement"]
        #total_template_calls += len(c.template_calls)
        #total_complement_calls += len(c.complement_calls)

    for w in xrange(workers):
        p = Process(target=run_methyl_caller, args=(work_queue, done_queue))
        p.start()
        jobs.append(p)
        work_queue.put('STOP')

    for p in jobs:
        p.join()

    done_queue.put('STOP')

if __name__ == "__main__":
    sys.exit(main(sys.argv))
