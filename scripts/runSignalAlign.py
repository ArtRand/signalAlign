#!/usr/bin/env python
"""Run signal-to-reference alignments
"""
from __future__ import print_function
import sys
sys.path.append("../")
from signalAlignLib import *
from multiprocessing import Process, Queue, current_process, Manager
from serviceCourse.file_handlers import FolderHandler
from argparse import ArgumentParser
from random import shuffle


def parse_args():
    parser = ArgumentParser(description=__doc__)

    parser.add_argument('--file_directory', '-d', action='store',
                        dest='files_dir', required=False, type=str, default=None,
                        help="directory with fast5s for alignment")
    parser.add_argument('--ref', '-r', action='store',
                        dest='ref', required=True, type=str, help="reference sequence to align to, in FASTA")

    parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for template events, if you don't want the default")
    parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for complement events, if you don't want the default")
    parser.add_argument('--in_complement_hmm_2', '-C2', action='store', dest='in_C_Hmm2',
                        required=False, type=str, default=None,
                        help="input HMM for complement events, if you don't want the default")
    parser.add_argument('--templateHDP', '-tH', action='store', dest='templateHDP', default=None)
    parser.add_argument('--complementHDP', '-cH', action='store', dest='complementHDP', default=None)

    parser.add_argument('--stateMachineType', '-smt', action='store', dest='stateMachineType', type=str,
                        required=True, default="threeState", help="decide which model to use, threeState by default")
    parser.add_argument('--threshold', '-t', action='store', dest='threshold', type=float, required=False,
                        default=None, help="posterior match probability threshold")
    parser.add_argument('--diagonalExpansion', '-e', action='store', dest='diag_expansion', type=int,
                        required=False, default=None, help="number of diagonals to expand around each anchor")
    parser.add_argument('--constraintTrim', '-m', action='store', dest='constraint_trim', type=int,
                        required=False, default=None, help='amount to remove from an anchor constraint')
    parser.add_argument('--target_regions', '-q', action='store', dest='target_regions', type=str,
                        required=False, default=None, help="tab separated table with regions to align to")
    parser.add_argument('---un-banded', '-ub', action='store_false', dest='banded',
                        default=True, help='flag, turn off banding')
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=4, type=int, help="number of jobs to run concurrently")
    parser.add_argument('-nb_files', '-n', action='store', dest='nb_files', required=False,
                        default=50, type=int, help="maximum number of reads to align")
    parser.add_argument('--cytosine_substitution', '-cs', action='append', default=None,
                        dest='cytosine_sub', required=False, type=str,
                        help="mutate cytosines to this letter in the reference")

    parser.add_argument('--sparse_output', '-s', action='store_true', default=False, dest='sparse',
                        help="Sparse output flag")
    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the alignments")

    args = parser.parse_args()
    return args


def aligner(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            alignment = SignalAlignment(**f)
            alignment.run()
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def main(args):
    # parse args
    args = parse_args()

    start_message = """
#   Starting Signal Align
#   Aligning files from: {fileDir}
#   Aligning to reference: {reference}
#   Aligning {nbFiles}
#   Using model: {model}
#   Using banding: {banding}
#   Aligning to regions in: {regions}
#   Input template HMM: {inThmm}
#   Input complement HMM: {inChmm}
    """.format(fileDir=args.files_dir, reference=args.ref, nbFiles=args.nb_files, banding=args.banded,
               inThmm=args.in_T_Hmm, inChmm=args.in_C_Hmm, model=args.stateMachineType, regions=args.target_regions)

    print(start_message, file=sys.stdout)

    if not os.path.isfile(args.ref):
        print("Did not find valid reference file", file=sys.stderr)
        sys.exit(1)

    # make directory to put temporary files
    temp_folder = FolderHandler()
    temp_dir_path = temp_folder.open_folder(args.out + "tempFiles_alignment")
    reference_seq = temp_folder.add_file_path("reference_seq.txt")
    make_temp_sequence(args.ref, True, reference_seq)

    # index the reference for bwa
    print("signalAlign - indexing reference", file=sys.stderr)
    bwa_ref_index = get_bwa_index(args.ref, temp_dir_path)
    print("signalAlign - indexing reference, done", file=sys.stderr)

    # parse the target regions, if provided
    if args.target_regions is not None:
        target_regions = TargetRegions(args.target_regions)
    else:
        target_regions = None

    workers = args.nb_jobs
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()
    jobs = []

    fast5s = [x for x in os.listdir(args.files_dir) if x.endswith(".fast5")]

    cytosine_sub = args.cytosine_sub[0] if args.cytosine_sub is not None else None

    nb_files = args.nb_files
    if nb_files < len(fast5s):
        shuffle(fast5s)
        fast5s = fast5s[:nb_files]

    for fast5 in fast5s:
        alignment_args = {
            "reference": reference_seq,
            "cytosine_substitution": cytosine_sub,
            "destination": temp_dir_path,
            "stateMachineType": args.stateMachineType,
            "bwa_index": bwa_ref_index,
            "in_templateHmm": args.in_T_Hmm,
            "in_complementHmm": args.in_C_Hmm,
            "in_complementHmm_pop1": args.in_C_Hmm2,
            "in_templateHdp": args.templateHDP,
            "in_complementHdp": args.complementHDP,
            "banded": args.banded,
            "sparse_output": args.sparse,
            "in_fast5": args.files_dir + fast5,
            "threshold": args.threshold,
            "diagonal_expansion": args.diag_expansion,
            "constraint_trim": args.constraint_trim,
            "target_regions": target_regions,
        }
        #alignment = SignalAlignment(**alignment_args)
        #alignment.run()
        work_queue.put(alignment_args)

    for w in xrange(workers):
        p = Process(target=aligner, args=(work_queue, done_queue))
        p.start()
        jobs.append(p)
        work_queue.put('STOP')

    for p in jobs:
        p.join()

    done_queue.put('STOP')
    print("\n#  signalAlign - finished alignments\n", file=sys.stderr)
    print("\n#  signalAlign - finished alignments\n", file=sys.stdout)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
