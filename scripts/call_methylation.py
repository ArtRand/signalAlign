#!/usr/bin/env python
"""Takes alignments from signalAlign and calls methylation status"""
from __future__ import print_function, division
import sys
sys.path.append("../")
from argparse import ArgumentParser
from alignmentAnalysisLib import CallMethylation
from signalAlignLib import parse_substitution_file, degenerate_enum
from variantCallingLib import get_alignments_labels_and_mask, get_reference_sequence
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
    parser.add_argument('--error_correct', action='store', default=None, required=False, type=int,
                        dest='error_correct', help="Enable error correction, provide error correction position")
    parser.add_argument('--degenerate', '-x', action='store', dest='degenerate', default="variant",
                        help="Specify degenerate nucleotide options: "
                             "variant -> {ACGT}, twoWay -> {CE} threeWay -> {CEO}")
    parser.add_argument('-n', required=False, action='store', type=int, dest='n', default=100,
                        help='Max number of alignments from each category to look at')
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=4, type=int, help="number of jobs to run concurrently")
    parser.add_argument('--out', '-o', action='store', type=str, required=True, dest='out')

    return parser.parse_args()


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

    if args.positions is not None and args.error_correct is not None:
        positions = {}
        f, b = parse_substitution_file(args.positions)
        positions['forward'] = f[1]
        positions['backward'] = b[1]
    elif args.error_correct is not None:
        positions = {'forward': range(args.error_correct, len(reference_sequence)),
                     'backward': range(args.error_correct, len(reference_sequence))
                     }
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
            "out_file": out_file,
            "positions": positions,
            "degenerate_type": degenerate_enum(args.degenerate),
        }
        #c = CallMethylation(**call_methyl_args)
        #c.write()
        work_queue.put(call_methyl_args)

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
