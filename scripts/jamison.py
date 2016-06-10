#!/usr/bin/env python
"""Run signal-to-reference alignments
"""
from __future__ import print_function
from signalAlignLib import *
from variantCallingLib import scan_for_proposals, update_reference_with_marginal_probs
from serviceCourse.file_handlers import FolderHandler
from argparse import ArgumentParser
from random import shuffle
from shutil import copyfile
from operator import itemgetter

STEP = 6


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
                             "variant -> {ACGT}, twoWay -> {CE} threeWay -> {CEO}")
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
                        default=4, type=int, help="number of jobs to run concurrently")
    parser.add_argument('--nb_files', '-n', action='store', dest='nb_files', required=False,
                        default=500, type=int, help="maximum number of reads to align")

    parser.add_argument('--cycles', dest='cycles', default=1, required=False, type=int)  # todo help string
    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the alignments")

    parser.add_argument('--corrected', dest='corrected', required=False, default='corrected.fa')  # todo help string

    args = parser.parse_args()
    return args


def group_sites_in_window(sites, window=6):
    def collect_group(start):
        i = start
        g = [sites[start]]
        while sites[i + 1][0] - sites[i][0] < window:
            g.append(sites[i + 1])
            i += 1
            if len(sites) <= i + 1:
                break
        return g, i + 1

    sites.sort()
    groups = []
    i = 0
    while i + 1 < len(sites):
        g, i = collect_group(i)
        groups.append(g)
    return [max(x, key=itemgetter(1))[0] for x in groups]


def get_first_sequence(input_fasta):
    input_sequence = ""
    for header, comment, sequence in read_fasta(input_fasta):
        input_sequence += sequence
        break
    return input_sequence


def main(args):
    # parse args
    args = parse_args()

    command_line = " ".join(sys.argv[:])
    print("Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)

    start_message = """
#   Starting Jamison Error-Correction
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
#   Performing {cycles} cycles
    """.format(fileDir=args.files_dir, reference=args.ref, nbFiles=args.nb_files, banding=args.banded,
               inThmm=args.in_T_Hmm, inChmm=args.in_C_Hmm, model=args.stateMachineType, regions=args.target_regions,
               tHdp=args.templateHDP, cHdp=args.complementHDP, cycles=args.cycles)

    print(start_message, file=sys.stdout)

    if not os.path.isfile(args.ref):
        print("Did not find valid reference file", file=sys.stderr)
        sys.exit(1)

    temp_folder = FolderHandler()
    temp_dir_path = temp_folder.open_folder(args.out + "tempFiles_errorCorrection")

    # initialize to input fasta
    reference_sequence_path = args.ref

    # list of alignment files
    fast5s = cull_fast5_files(args.files_dir, args.nb_files)

    for cycle in range(0, args.cycles):
        # index the reference for bwa this is a string with the path to the index
        bwa_ref_index = get_bwa_index(reference_sequence_path, temp_dir_path)

        # unpack the reference sequence
        reference_sequence_string = get_first_sequence(reference_sequence_path)

        alignment_args = {
            "path_to_EC_refs": None,
            "destination": temp_dir_path,
            "stateMachineType": args.stateMachineType,
            "bwa_index": bwa_ref_index,
            "in_templateHmm": args.in_T_Hmm,
            "in_complementHmm": args.in_C_Hmm,
            "in_templateHdp": args.templateHDP,
            "in_complementHdp": args.complementHDP,
            "banded": args.banded,
            "sparse_output": True,
            "threshold": args.threshold,
            "diagonal_expansion": args.diag_expansion,
            "constraint_trim": args.constraint_trim,
            "target_regions": None,
            "degenerate": degenerate_enum(args.degenerate),
        }

        proposals = scan_for_proposals(temp_folder, STEP, reference_sequence_string, fast5s, alignment_args,
                                       args.nb_jobs)

        proposals = group_sites_in_window(proposals, 6)

        print("Cycle {cycle} - Got {nb} sites to check: {sites}".format(nb=len(proposals),
                                                                        sites=proposals,
                                                                        cycle=cycle))

        updated_reference_string = update_reference_with_marginal_probs(temp_folder, proposals,
                                                                        reference_sequence_string, fast5s,
                                                                        alignment_args, args.nb_jobs)

        updated_reference_path = temp_folder.add_file_path("cycle_snapshot.{cycle}.fa".format(cycle=cycle))

        write_fasta("jamison{}".format(cycle), updated_reference_string, open(updated_reference_path, 'w'))

        reference_sequence_path = updated_reference_path

    # copy final file
    copyfile(reference_sequence_path, temp_dir_path + args.corrected)

    return

if __name__ == "__main__":
    sys.exit(main(sys.argv))
