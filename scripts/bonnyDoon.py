#!/usr/bin/env python
"""Run signal-to-reference alignments
"""
from __future__ import print_function
import pandas as pd

from signalAlignLib import *
from alignmentAnalysisLib import CallMethylation, get_first_sequence
from variantCallingLib import scan_for_proposals
from multiprocessing import Process, Queue, current_process, Manager
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
    # todo help string
    parser.add_argument('--cycles', dest='cycles', default=1, required=False, type=int)

    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the alignments")
    # todo help string
    parser.add_argument('--corrected', dest='corrected', required=False, default='corrected.fa')

    args = parser.parse_args()
    return args


def group_sites_in_window2(sites, window=6):
    def collect_group(start):
        i = start
        g = [sites[start]]
        while sites[i + 1] - sites[i] < window:
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
    return groups


def make_degenerate_reference(input_fasta, start, forward_sequence_path, backward_sequence_path,
                              block_size=1, step=6):
    """
    input_sequence: string, input nucleotide sequence
    out_path: string, path to directory to put new sequences with substituted degenerate characters
    block_size: not implemented
    step: number of bases between degenerate characters
    :return (subbed sequence, complement subbed sequence)
    """

    input_sequence = get_first_sequence(input_fasta)

    complement_sequence = reverse_complement(dna=input_sequence, reverse=False, complement=True)

    t_seq = list(input_sequence)
    c_seq = list(complement_sequence)

    positions = xrange(start, len(input_sequence), step)
    for position in positions:
        t_seq[position] = "X"
        c_seq[position] = "X"

    t_seq = ''.join(t_seq)
    c_seq = ''.join(c_seq)

    sequence_length = len(t_seq)

    with open(forward_sequence_path, 'w') as f:
        f.write("{seq}".format(seq=t_seq))
    with open(backward_sequence_path, 'w') as f:
        f.write("{seq}".format(seq=c_seq))

    return True, sequence_length


def aligner(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            alignment = SignalAlignment(**f)
            alignment.run()
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def run_methyl_caller(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            c = CallMethylation(**f)
            c.write()
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def load_data(file_path):
    data = pd.read_table(file_path,
                         usecols=(0, 1, 2, 3, 4, 5, 6),
                         names=['site', 'strand', 'pA', 'pC', 'pG', 'pT', 'read'],
                         dtype={'site': np.int64,
                                'strand': np.str,
                                'pC': np.float64,
                                'pmC': np.float64,
                                'phmC': np.float64,
                                'read': np.str,
                                })
    return data


def symbol_to_base(symbol):
    return ["A", "C", "G", "T"][symbol]


def rc_probs(probs):
    return [probs[3], probs[2], probs[1], probs[0]]


def update_reference(data, reference_sequence, min_depth=0, get_sites=False):
    d = load_data(data)

    ref = get_first_sequence(reference_sequence)
    ref = list(ref)

    candidate_sites = []
    add_to_candidates = candidate_sites.append

    for g, x in d.groupby("site"):
        marginal_forward_p = pd.Series(0, ['pA', 'pC', 'pG', 'pT'])
        marginal_backward_p = pd.Series(0, ['pA', 'pC', 'pG', 'pT'])
        assert(len(x['site'].unique()) == 1)
        site = x['site'].unique()[0]

        if len(x['read']) < min_depth:
            continue

        for i, read in x.iterrows():
            if ((read['read'].endswith(".forward.tsv") and read['strand'] == 't') or
                    (read['read'].endswith(".backward.tsv") and read['strand'] == 'c')):
                direction = True
            else:
                direction = False

            if direction:
                marginal_forward_p += read[['pA', 'pC', 'pG', 'pT']]
            else:
                marginal_backward_p += read[['pA', 'pC', 'pG', 'pT']]

        marginal_prob = marginal_forward_p + rc_probs(marginal_backward_p)

        normed_marginal_probs = marginal_prob.map(lambda x: x / sum(marginal_prob))
        called_base = normed_marginal_probs.argmax()[1]
        #called_base = marginal_prob.map(lambda x: x / sum(marginal_prob)).argmax()[1]

        if called_base != ref[site]:
            if get_sites is False:
                print("Changing {orig} to {new} at {site}".format(orig=ref[site], new=called_base, site=site))
                ref[site] = called_base
            else:
                print("Proposing edit at {site} from {orig} to {new}, \n{probs}"
                      "".format(orig=ref[site], new=called_base, site=site, probs=normed_marginal_probs))
                difference = normed_marginal_probs.max() - normed_marginal_probs["p" + ref[site]]
                print(difference)
                add_to_candidates((site, difference))

    if get_sites is True:
        return candidate_sites
    else:
        return ''.join(ref)


def main(args):
    # parse args
    args = parse_args()

    command_line = " ".join(sys.argv[:])
    print("Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)

    start_message = """
#   Starting BonnyDoon Error-Correction
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
    # cull the MinION files
    fast5s = cull_fast5_files(args.files_dir, args.nb_files)

    # get the (input) reference sequence
    if not os.path.isfile(args.ref):
        print("Did not find valid reference file", file=sys.stderr)
        sys.exit(1)
    reference_sequence_path = args.ref

    # unpack the reference sequence
    reference_sequence_string = get_first_sequence(reference_sequence_path)

    # make a working folder in the specified directory
    temp_folder = FolderHandler()
    temp_dir_path = temp_folder.open_folder(args.out + "tempFiles_errorCorrection")

    # index the reference for bwa this is a string with the path to the index
    bwa_ref_index = get_bwa_index(reference_sequence_path, temp_dir_path)

    # alignment args are the parameters to the HMM/HDP model, and don't change
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

    # get the sites that have proposed edits
    proposals = scan_for_proposals(temp_folder, STEP, reference_sequence_string, fast5s, alignment_args, args.nb_jobs)
    proposals = group_sites_in_window2([x[0] for x in proposals], 6)






    return

if __name__ == "__main__":
    sys.exit(main(sys.argv))
