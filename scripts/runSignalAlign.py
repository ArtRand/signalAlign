#!/usr/bin/env python
"""Main driver script for running an ionic current-to-sequence alignment on a single machine.
"""
from __future__ import print_function

import sys
import os

#import pysam

from argparse import ArgumentParser
from random import shuffle
from multiprocessing import Process, current_process, Manager

from signalalign.signalAlignment import SignalAlignment
from signalalign.utils import processReferenceFasta, parseFofn
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.utils.bwaWrapper import getBwaIndex
from signalalign.motif import getDegenerateEnum


def signalAlignSourceDir():
    return "/".join(os.path.abspath(__file__).split("/")[:-1])  # returns path without runSignalAlign


def resolvePath(p):
    if p is None:
        return None
    elif p.startswith("/"):
        return p
    else:
        return os.path.abspath(p)


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--file_directory', '-d', action='store',
                        dest='files_dir', required=True, type=str, default=None,
                        help="directory with MinION fast5 reads to align")
    parser.add_argument('--ref', '-r', action='store',
                        dest='ref', required=True, type=str,
                        help="reference sequence to align to, in FASTA")
    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the alignments")
    # optional arguments
    parser.add_argument("--2d", action='store_true', dest="twoD", default=False, help="flag, specify if using 2D reads")
    parser.add_argument("--bwt", action='store', dest="bwt", default=None, help="path to BWT files. example: ../ref.fasta")
    parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for template events, if you don't want the default")
    parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for complement events, if you don't want the default")
    parser.add_argument('--template_hdp', '-tH', action='store', dest='templateHDP', default=None,
                        help="template serialized HDP file")
    parser.add_argument('--complement_hdp', '-cH', action='store', dest='complementHDP', default=None,
                        help="complement serialized HDP file")
    parser.add_argument('--degenerate', '-x', action='store', dest='degenerate', default="variant",
                        help="Specify degenerate nucleotide options: "
                             "variant -> {ACGT}, cytosine2 -> {CE} cytosine3 -> {CEO} adenosine -> {AI}")
    parser.add_argument('--stateMachineType', '-smt', action='store', dest='stateMachineType', type=str,
                        default="threeState", help="decide which model to use, threeState by default")
    parser.add_argument('--file_of_files', '-fofn', action='store', dest='fofn', required=False, type=str, default=None,
                        help="text file containing absolute paths to files to use")
    parser.add_argument('--threshold', '-t', action='store', dest='threshold', type=float, required=False,
                        default=0.01, help="posterior match probability threshold, Default: 0.01")
    parser.add_argument('--diagonalExpansion', '-e', action='store', dest='diag_expansion', type=int,
                        required=False, default=None,
                        help="number of diagonals to expand around each anchor default: 50")
    parser.add_argument('--constraintTrim', '-m', action='store', dest='constraint_trim', type=int,
                        required=False, default=None, help='amount to remove from an anchor constraint')
    parser.add_argument('--target_regions', '-q', action='store', dest='target_regions', type=str,
                        required=False, default=None, help="tab separated table with regions to align to")
    parser.add_argument("--motif", action="store", dest="motif_key", default=None)
    parser.add_argument('--ambiguity_positions', '-p', action='store', required=False, default=None,
                        dest='ambiguity_positions', help="Ambiguity positions")
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=4, type=int, help="number of jobs to run in parallel")
    parser.add_argument('--nb_files', '-n', action='store', dest='nb_files', required=False,
                        default=500, type=int, help="maximum number of reads to align")
    parser.add_argument('--ambig_char', '-X', action='store', required=False, default=None, type=str, dest='ambig_char',
                        help="")
    parser.add_argument('--output_format', '-f', action='store', default="full", dest='outFmt',
                        help="output format: full, variantCaller, or assignments. Default: full")
    parser.add_argument('--debug', action='store_true', dest="DEBUG", default=False)

    args = parser.parse_args()
    return args


def aligner(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            alignment = SignalAlignment(**f)
            alignment.run()
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def concat_variant_call_files(path):
    concat_command = "cat {path}/*.tsv > {path}/probs.tsv".format(path=path)
    os.system(concat_command)
    return


def main(args):
    # parse args
    args = parse_args()

    command_line = " ".join(sys.argv[:])
    print("Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)

    # get absolute paths to inputs
    args.files_dir           = resolvePath(args.files_dir)
    args.ref                 = resolvePath(args.ref)
    args.out                 = resolvePath(args.out)
    args.bwt                 = resolvePath(args.bwt)
    args.in_T_Hmm            = resolvePath(args.in_T_Hmm)
    args.in_C_Hmm            = resolvePath(args.in_C_Hmm)
    args.templateHDP         = resolvePath(args.templateHDP)
    args.complementHDP       = resolvePath(args.complementHDP)
    args.fofn                = resolvePath(args.fofn)
    args.target_regions      = resolvePath(args.target_regions)
    args.ambiguity_positions = resolvePath(args.ambiguity_positions)
    start_message = """
#   Starting Signal Align
#   Aligning files from: {fileDir}
#   Aligning to reference: {reference}
#   Aligning maximum of {nbFiles} files
#   Using model: {model}
#   Using banding: True
#   Aligning to regions in: {regions}
#   Non-default template HMM: {inThmm}
#   Non-default complement HMM: {inChmm}
#   Template HDP: {tHdp}
#   Complement HDP: {cHdp}
    """.format(fileDir=args.files_dir, reference=args.ref, nbFiles=args.nb_files,
               inThmm=args.in_T_Hmm, inChmm=args.in_C_Hmm, model=args.stateMachineType, regions=args.target_regions,
               tHdp=args.templateHDP, cHdp=args.complementHDP)

    print(start_message, file=sys.stdout)

    if args.files_dir is None and args.fofn is None:
        print("Need to provide directory with .fast5 files of fofn", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(args.ref):
        print("Did not find valid reference file, looked for it {here}".format(here=args.ref), file=sys.stderr)
        sys.exit(1)

    # make directory to put temporary files
    temp_folder = FolderHandler()
    temp_dir_path = temp_folder.open_folder(args.out + "/tempFiles_alignment")

    reference_map = processReferenceFasta(fasta=args.ref,
                                          motif_key=args.motif_key,
                                          work_folder=temp_folder,
                                          sub_char=args.ambig_char,
                                          positions_file=args.ambiguity_positions)

    # index the reference for bwa
    if args.bwt is not None:
        print("[RunSignalAlign]NOTICE - using provided BWT %s" % args.bwt)
        bwa_ref_index = args.bwt
    else:
        print("signalAlign - indexing reference", file=sys.stderr)
        bwa_ref_index = getBwaIndex(args.ref, temp_dir_path)
        print("signalAlign - indexing reference, done", file=sys.stderr)

    # setup workers for multiprocessing
    workers = args.nb_jobs
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()
    jobs = []

    # list of read files
    if args.fofn is not None:
        fast5s = [x for x in parseFofn(args.fofn) if x.endswith(".fast5")]
    else:
        fast5s = ["/".join([args.files_dir, x]) for x in os.listdir(args.files_dir) if x.endswith(".fast5")]

    nb_files = args.nb_files
    if nb_files < len(fast5s):
        shuffle(fast5s)
        fast5s = fast5s[:nb_files]

    # change paths to the source directory
    os.chdir(signalAlignSourceDir())

    print("[runSignalAlign]:NOTICE: Got {} files to align".format(len(fast5s)), file=sys.stdout)
    for fast5 in fast5s:
        print(fast5)
        alignment_args = {
            "reference_map": reference_map,
            "destination": temp_dir_path,
            "stateMachineType": args.stateMachineType,
            "bwa_index": bwa_ref_index,
            "in_templateHmm": args.in_T_Hmm,
            "in_complementHmm": args.in_C_Hmm,
            "in_templateHdp": args.templateHDP,
            "in_complementHdp": args.complementHDP,
            "output_format": args.outFmt,
            "in_fast5": fast5,
            "threshold": args.threshold,
            "diagonal_expansion": args.diag_expansion,
            "constraint_trim": args.constraint_trim,
            "degenerate": getDegenerateEnum(args.degenerate),
            "twoD_chemistry": args.twoD,
            "target_regions": args.target_regions,
        }
        if args.DEBUG:
            alignment = SignalAlignment(**alignment_args)
            alignment.run()
        else:
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
