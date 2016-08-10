#!/usr/bin/env python
"""Master pipeline script for generating trained HDPs for MinION signal data

Input: alignments using non-HDP model
Output: trained HDP and model

The objective of this pipeline is to:
    1. use input alignments to make 'build alignment' for generating the initial HDP
    2. generates the initial HDP
    3. trains the HDP on MinION reads
    4. outputs distributions for all kmers from the HDP
"""

import os
import sys
from argparse import ArgumentParser
from subprocess import check_call, Popen
from shutil import copyfile


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # build alignment
    parser.add_argument('--C_alignments', '-C', action='store',
                        dest='C_alns', required=False, type=str, default=None,
                        help="C files")
    parser.add_argument('--mC_alignments', '-mC', action='store',
                        dest='mC_alns', required=False, type=str, default=None,
                        help="mC files")
    parser.add_argument('--hmC_alignments', '-hmC', action='store',
                        dest='hmC_alns', required=False, type=str, default=None,
                        help="hmC files")
    parser.add_argument('--number_of_assignments', '-n', action='store', type=int, default=10000,
                        dest='max_assignments',
                        help='total number of assignments to collect FOR EACH GROUP')
    # initial HDP
    parser.add_argument('--build_alignment', action='store', type=str, default=None,
                        required=False, dest='build_alignment')
    parser.add_argument('--threshold', '-t', action='store', type=float, default=0.0, dest='threshold')
    parser.add_argument('--hdp_type', action='store', type=str, required=False, dest='hdp_type', default='Prior',
                        help="Build Hdp, specify type, options: "
                             "Prior, Fixed, twoWay. twoWay is a Prior-type model (recommended)")
    parser.add_argument('--template_model', '-tM', action='store', type=str, dest='template_lookup',
                        required=True, help="Input template lookup table")
    parser.add_argument('--complement_model', '-cM', action='store', type=str, dest='complement_lookup',
                        required=True, help="Input complement lookup table")
    # fixed concentration models
    parser.add_argument('--base_gamma', '-B', action='store', type=float, default=1.0, dest='base_gamma',
                        required=False)
    parser.add_argument('--middle_gamma', '-M', action='store', type=float, default=1.0, dest='middle_gamma',
                        required=False)
    parser.add_argument('--leaf_gamma', '-L', action='store', type=float, default=1.0, dest='leaf_gamma',
                        required=False)
    # gamma prior models
    parser.add_argument('--base_alpha', '-Ba', action='store', type=float, default=1.0, dest='base_alpha',
                        required=False)
    parser.add_argument('--base_beta', '-Bb', action='store', type=float, default=1.0, dest='base_beta',
                        required=False)
    parser.add_argument('--middle_alpha', '-Ma', action='store', type=float, default=1.0, dest='middle_alpha',
                        required=False)
    parser.add_argument('--middle_beta', '-Mb', action='store', type=float, default=1.0, dest='middle_beta',
                        required=False)
    parser.add_argument('--leaf_alpha', '-La', action='store', type=float, default=1.0, dest='leaf_alpha',
                        required=False)
    parser.add_argument('--leaf_beta', '-Lb', action='store', type=float, default=1.0, dest='leaf_beta',
                        required=False)
    # gibbs
    parser.add_argument('--samples', '-s', action='store', type=int, default=10000, dest='gibbs_samples')
    parser.add_argument('--thinning', '-th', action='store', type=int, default=100, dest='thinning')
    parser.add_argument('--verbose', action='store_true', default=False, dest='verbose')
    # sample grid
    parser.add_argument('--grid_start', action='store', type=float, default=30.0, dest='grid_start')
    parser.add_argument('--grid_end', action='store', type=float, default=90.0, dest='grid_end')
    parser.add_argument('--grid_length', action='store', type=int, default=1200, dest='grid_length')

    parser.add_argument('--out', '-o', action='store', type=str, required=True, dest='out')

    return parser.parse_args()


def get_set_of_hdp_types(request):
    if request == 'prior':
        return [1, 3, 5, 7, 9]
    else:
        return [1, 2, 4, 6, 8]


def get_hdp_type(requested_type):
        hdp_types = {
            "singleLevelFixed": 0,
            "singleLevelPrior": 1,
            "multisetFixed": 2,
            "multisetPrior": 3,
            "compFixed": 4,
            "compPrior": 5,
            "middleNtsFixed": 6,
            "middleNtsPrior": 7,
            "groupMultisetFixed": 8,
            "groupMultisetPrior": 9,
            "singleLevelPrior2": 10,
            "multisetPrior2": 11,
            "multisetPriorEcoli": 12,
            "singleLevelPriorEcoli": 13
        }
        assert (requested_type in hdp_types.keys()), "Requested HDP type is invalid, got {}".format(requested_type)
        return hdp_types[requested_type]


def count_lines_in_build_alignment(build_alignment_path):
    count = 0
    for line in open(build_alignment_path, 'r').xreadlines():
        count += 1
    return count


def kmer_length_from_model(template_model_file, complement_model_file):
    def get_kmer_length(model_file):
        with open(model_file, "r") as fH:
            line = fH.readline().split()
            assert len(line) == 4, "HdpPipeline ERROR: wrong header in model file {}".format(model_file)
            kmer_length = int(line[3])
            fH.close()
        return kmer_length
    template_kmer_length = get_kmer_length(template_model_file)
    complement_kmer_length = get_kmer_length(complement_model_file)
    assert template_kmer_length == complement_kmer_length
    return template_kmer_length


def get_initial_hdp_args(args, hdp_type):
    # if we're making a HDP with fixed concentration parameters
    if hdp_type in [0, 2, 4, 6, 8]:
        assert None not in [args.base_gamma, args.leaf_gamma], \
            "ERROR: need to specify concentration parameters for type {}".format(hdp_type)
        if hdp_type == 0:
            return "-B {base} -L {leaf} ".format(base=args.base_gamma, leaf=args.leaf_gamma)
        else:
            assert args.middle_gamma is not None, "ERROR: need to specify middle concentration param"
            return "-B {base} -M {middle} -L {leaf} ".format(base=args.base_gamma, middle=args.middle_gamma,
                                                             leaf=args.leaf_gamma)
    else:
        assert None not in [args.base_alpha, args.base_beta, args.leaf_alpha, args.leaf_beta], \
                "ERROR: missing Gamma prior hyper parameters"
        if hdp_type == 1 or hdp_type == 10:
            return "-g {Ba} -r {Bb} -i {La} -u {Lb} ".format(Ba=args.base_alpha, Bb=args.base_beta,
                                                             La=args.leaf_alpha, Lb=args.leaf_beta)
        else:
            assert None not in [args.middle_alpha, args.middle_beta], "ERROR: need middle hyper parameters"
            return "-g {Ba} -r {Bb} -j {Ma} -y {Mb} -i {La} -u {Lb} ".format(Ba=args.base_alpha, Bb=args.base_beta,
                                                                             Ma=args.middle_alpha, Mb=args.middle_beta,
                                                                             La=args.leaf_alpha, Lb=args.leaf_beta)


# globals
HDP_TYPES = [
    ("singleLevelFixed", 0),
    ("singleLevelPrior", 1),
    ("multisetFixed", 2),
    ("multisetPrior", 3),
    ("compFixed", 4),
    ("compPrior", 5),
    ("middleNtsFixed", 6),
    ("middleNtsPrior", 7),
    ("groupMultisetFixed", 8),
    ("groupMultisetPrior", 9),
]

HDP_TYPES_2 = [
    ("singleLevelPrior2", 10),
    ("multisetPrior2", 11),
]

HDP_TYPES_ECOLI = [
    ("multisetPriorEcoli", 12),
    ("singleLevelPriorEcoli", 13),
]

# Pipeline Script
args = parse_args()  # parse arguments
working_directory = args.out  # this is the directory we will use for everything
assert os.path.isdir(working_directory), "ERROR: the working directory you specified doesn't exist."
pipeline_log = open(working_directory + "pipeline.log", 'a')
command_line = " ".join(sys.argv[:])
pipeline_log.write("[pipeline] Command Line: {}\n".format(command_line))
signalAlign_directory = "../../signalAlign/"
build_alignment_location = working_directory + "buildAlignment.tsv"
if args.build_alignment is None:
    # build alignment
    build_alignment_command = "{sA}scripts/makeBuildAlignments.py -o={bA} -t={threshold} -n={nbAssignments} " \
                              "".format(sA=signalAlign_directory, C=args.C_alns, mC=args.mC_alns,
                                        threshold=args.threshold, hmC=args.hmC_alns, bA=build_alignment_location,
                                        nbAssignments=args.max_assignments)
    approx_total_build_assignments = 0  # keep track of about how many assignments we're going to get for gibbs burn in
    if args.C_alns is not None:  # add the alignments to the command
        build_alignment_command += "-C={C} ".format(C=args.C_alns)
        approx_total_build_assignments += args.max_assignments
    if args.mC_alns is not None:
        build_alignment_command += "-mC={mC} ".format(mC=args.mC_alns)
        approx_total_build_assignments += args.max_assignments
    if args.hmC_alns is not None:
        build_alignment_command += "-hmC={hmC} ".format(hmC=args.hmC_alns)
        approx_total_build_assignments += args.max_assignments
    pipeline_log.write("[pipeline] NOTICE: Making build alignment using files from:\n\t{C}\n\t{mC}\n\t{hmC}\n"
                       "".format(C=args.C_alns, mC=args.mC_alns, hmC=args.hmC_alns))
    pipeline_log.write("[pipeline] Command: {}\n".format(build_alignment_command))
    check_call(build_alignment_command.split(), stderr=pipeline_log, stdout=pipeline_log)
else:
    pipeline_log.write("[pipeline] NOTICE: using build alignment {}".format(args.build_alignment))
    assert os.path.isfile(args.build_alignment), "ERROR: Didn't find input BuildAlignment"
    copyfile(args.build_alignment, build_alignment_location)
    approx_total_build_assignments = count_lines_in_build_alignment(build_alignment_location)

# initial HDP
assert (os.path.isfile(build_alignment_location)), "ERROR: Didn't find build alignment"
assert (os.path.exists("./buildHdpUtil")), "ERROR: Didn't find buildHdpUtil"
pipeline_log.write("[pipeline] NOTICE: Making initial HDP of type {}\n".format(args.hdp_type))

initial_hdp_build_out = open(working_directory + "build_initial_hdp.out", 'w')
initial_hdp_build_err = open(working_directory + "build_initial_hdp.err", 'w')
kmer_length = kmer_length_from_model(args.template_lookup, args.complement_lookup)
template_lookup_table = " -T" + args.template_lookup
complement_lookup_table = " -C" + args.complement_lookup
verbose_flag = "--verbose " if args.verbose is True else ""
build_commands = []
if args.hdp_type == "cytosine2":
    hdp_types = HDP_TYPES_2
elif args.hdp_type == "ecoli":
    hdp_types = HDP_TYPES_ECOLI
else:
    hdp_types = HDP_TYPES[1::2] if args.hdp_type == "Prior" else HDP_TYPES[::2]
for hdp_type, i, in hdp_types:
    template_hdp_location = working_directory + "template." + hdp_type + ".nhdp"
    complement_hdp_location = working_directory + "complement." + hdp_type + ".nhdp"
    build_initial_hdp_command = "./buildHdpUtil {verbose}-p {hdpType} -v {tHdpLoc} -w {cHdpLoc} -l {buildAln} " \
                                "-a {kmerLength} -n {samples} -I {burnIn} -t {thin} -s {start} -e {end} " \
                                "-k {len}{tL}{cL} " \
                                "".format(hdpType=i, tHdpLoc=template_hdp_location,
                                          cHdpLoc=complement_hdp_location, buildAln=build_alignment_location,
                                          samples=args.gibbs_samples, burnIn=32 * approx_total_build_assignments,
                                          thin=args.thinning, start=args.grid_start, end=args.grid_end,
                                          len=args.grid_length, verbose=verbose_flag, tL=template_lookup_table,
                                          cL=complement_lookup_table, kmerLength=kmer_length)
    build_initial_hdp_command += get_initial_hdp_args(args=args, hdp_type=i)
    build_commands.append(build_initial_hdp_command)
    pipeline_log.write("[pipeline] Command: {}\n".format(build_initial_hdp_command))

procs = [Popen(x.split(), stdout=initial_hdp_build_out, stderr=initial_hdp_build_err) for x in build_commands]
status = [p.wait() for p in procs]

initial_hdp_build_out.close()
initial_hdp_build_err.close()

pipeline_log.write("[pipeline] DONE.\n")
pipeline_log.close()





