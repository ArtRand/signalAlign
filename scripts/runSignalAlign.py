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
                        dest='files_dir', required=True, type=str, default=None,
                        help="directory with MinION fast5 reads to align")
    parser.add_argument('--ref', '-r', action='store',
                        dest='ref', required=True, type=str,
                        help="reference sequence to align to, in FASTA")
    parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for template events, if you don't want the default")
    parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for complement events, if you don't want the default")
    parser.add_argument('--templateHDP', '-tH', action='store', dest='templateHDP', default=None,
                        help="template serialized HDP file")
    parser.add_argument('--complementHDP', '-cH', action='store', dest='complementHDP', default=None,
                        help="complement serialized HDP file")
    parser.add_argument('--degenerate', '-x', action='store', dest='degenerate', default="variant",
                        help="Specify degenerate nucleotide options: "
                             "variant -> {ACGT}, cytosine2 -> {CE} cytosine3 -> {CEO} adenosine {AI}")
    parser.add_argument('--stateMachineType', '-smt', action='store', dest='stateMachineType', type=str,
                        default="threeState", help="decide which model to use, threeState by default")
    parser.add_argument('--threshold', '-t', action='store', dest='threshold', type=float, required=False,
                        default=None, help="posterior match probability threshold, Default: 0.01")
    parser.add_argument('--diagonalExpansion', '-e', action='store', dest='diag_expansion', type=int,
                        required=False, default=None, help="number of diagonals to expand around each anchor")
    parser.add_argument('--constraintTrim', '-m', action='store', dest='constraint_trim', type=int,
                        required=False, default=None, help='amount to remove from an anchor constraint')
    parser.add_argument('--target_regions', '-q', action='store', dest='target_regions', type=str,
                        required=False, default=None, help="tab separated table with regions to align to")
    parser.add_argument('---un-banded', '-ub', action='store_false', dest='banded',
                        default=True, help='flag, turn off banding')
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=4, type=int, help="number of jobs to run in parallel")
    parser.add_argument('--nb_files', '-n', action='store', dest='nb_files', required=False,
                        default=500, type=int, help="maximum number of reads to align")
    parser.add_argument('--ambiguity_positions', '-p', action='store', required=False, default=None,
                        dest='substitution_file', help="Ambiguity positions")
    parser.add_argument('--ambig_char', '-X', action='store', required=False, default="X", type=str, dest='ambig_char',
                        help="Character to substitute at positions, default is 'X'.")
    parser.add_argument('--output_format', '-f', action='store', default="full", dest='outFmt',
                        help="output format: full, variantCaller, or assignments. Default: full")
    parser.add_argument('--error_correct', action='store_true', default=False, required=False,
                        dest='error_correct', help="Enable error correction")
    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the alignments")

    args = parser.parse_args()
    return args


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

    command_line = " ".join(sys.argv[:])
    print("Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)

    start_message = """
#   Starting Signal Align
#   Aligning files from: {fileDir}
#   Aligning to reference: {reference}
#   Aligning maximum of {nbFiles} files
#   Using model: {model}
#   Using banding: {banding}
#   Aligning to regions in: {regions}
#   Non-default template HMM: {inThmm}
#   Non-default complement HMM: {inChmm}
#   Template HDP: {tHdp}
#   Complement HDP: {cHdp}
    """.format(fileDir=args.files_dir, reference=args.ref, nbFiles=args.nb_files, banding=args.banded,
               inThmm=args.in_T_Hmm, inChmm=args.in_C_Hmm, model=args.stateMachineType, regions=args.target_regions,
               tHdp=args.templateHDP, cHdp=args.complementHDP)

    print(start_message, file=sys.stdout)

    if not os.path.isfile(args.ref):
        print("Did not find valid reference file", file=sys.stderr)
        sys.exit(1)

    # make directory to put temporary files
    temp_folder = FolderHandler()
    temp_dir_path = temp_folder.open_folder(args.out + "tempFiles_alignment")

    if args.error_correct is True:
        write_degenerate_reference_set(input_fasta=args.ref, out_path=temp_dir_path)
        plus_strand_sequence = None
        minus_strand_sequence = None
    else:
        # parse the substitution file, if given
        plus_strand_sequence = temp_folder.add_file_path("forward_reference.txt")
        minus_strand_sequence = temp_folder.add_file_path("backward_reference.txt")
        if args.substitution_file is not None:
            add_ambiguity_chars_to_reference(input_fasta=args.ref,
                                             substitution_file=args.substitution_file,
                                             sequence_outfile=plus_strand_sequence,
                                             rc_sequence_outfile=minus_strand_sequence,
                                             degenerate_type=args.degenerate,
                                             ambig_char=args.ambig_char)
        else:
            make_temp_sequence(fasta=args.ref, sequence_outfile=plus_strand_sequence,
                               rc_sequence_outfile=minus_strand_sequence)

    # index the reference for bwa
    print("signalAlign - indexing reference", file=sys.stderr)
    bwa_ref_index = get_bwa_index(args.ref, temp_dir_path)
    print("signalAlign - indexing reference, done", file=sys.stderr)

    # parse the target regions, if provided
    # TODO make this the same as the 'labels' file
    if args.target_regions is not None:
        target_regions = TargetRegions(args.target_regions)
    else:
        target_regions = None

    # setup workers for multiprocessing
    workers = args.nb_jobs
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()
    jobs = []

    # list of read files
    fast5s = [x for x in os.listdir(args.files_dir) if x.endswith(".fast5")]

    nb_files = args.nb_files
    if nb_files < len(fast5s):
        shuffle(fast5s)
        fast5s = fast5s[:nb_files]
    print("Got {} files to align".format(len(fast5s)), file=sys.stdout)
    for fast5 in fast5s:
        alignment_args = {
            "forward_reference": plus_strand_sequence,
            "backward_reference": minus_strand_sequence,
            "path_to_EC_refs": (temp_dir_path if args.error_correct else None),
            "destination": temp_dir_path,
            "stateMachineType": args.stateMachineType,
            "bwa_index": bwa_ref_index,
            "in_templateHmm": args.in_T_Hmm,
            "in_complementHmm": args.in_C_Hmm,
            "in_templateHdp": args.templateHDP,
            "in_complementHdp": args.complementHDP,
            "banded": args.banded,
            "output_format": args.outFmt,
            "in_fast5": args.files_dir + fast5,
            "threshold": args.threshold,
            "diagonal_expansion": args.diag_expansion,
            "constraint_trim": args.constraint_trim,
            "target_regions": target_regions,
            "degenerate": degenerate_enum(args.degenerate),
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
