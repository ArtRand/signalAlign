#!/usr/bin/env python
"""Master pipeline script for generating trained HDPs for MinION signal data

Input: alignments using non-HDP model, desired HDP type
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
    parser.add_argument('--threshold', '-t', action='store', type=float, default=0.2, dest='threshold')
    parser.add_argument('--hdp_type', action='store', type=str, required=True, dest='hdp_type',
                        help="Build Hdp, specify type, options: "
                             "singleLevelFixed, singleLevelPrior, multisetFixed, multisetPrior")
    # fixed concentration models
    parser.add_argument('--base_gamma', '-B', action='store', type=float, default=None, dest='base_gamma',
                        required=False)
    parser.add_argument('--middle_gamma', '-M', action='store', type=float, default=None, dest='middle_gamma',
                        required=False)
    parser.add_argument('--leaf_gamma', '-L', action='store', type=float, default=None, dest='leaf_gamma',
                        required=False)
    # gamma prior models
    parser.add_argument('--base_alpha', '-Ba', action='store', type=float, default=None, dest='base_alpha',
                        required=False)
    parser.add_argument('--base_beta', '-Bb', action='store', type=float, default=None, dest='base_beta',
                        required=False)
    parser.add_argument('--middle_alpha', '-Ma', action='store', type=float, default=None, dest='middle_alpha',
                        required=False)
    parser.add_argument('--middle_beta', '-Mb', action='store', type=float, default=None, dest='middle_beta',
                        required=False)
    parser.add_argument('--leaf_alpha', '-La', action='store', type=float, default=None, dest='leaf_alpha',
                        required=False)
    parser.add_argument('--leaf_beta', '-Lb', action='store', type=float, default=None, dest='leaf_beta',
                        required=False)
    # gibbs
    parser.add_argument('--samples', '-s', action='store', type=int, default=10000, dest='gibbs_samples')
    parser.add_argument('--burnIn', '-I', action='store', type=int, default=100000, dest='burnIn')
    parser.add_argument('--thinning', '-th', action='store', type=int, default=100, dest='thinning')
    parser.add_argument('--verbose', action='store_true', default=False, dest='verbose')
    # sample grid
    parser.add_argument('--grid_start', action='store', type=float, default=30.0, dest='grid_start')
    parser.add_argument('--grid_end', action='store', type=float, default=90.0, dest='grid_end')
    parser.add_argument('--grid_length', action='store', type=int, default=1200, dest='grid_length')
    # train models
    parser.add_argument('--file_directory', '-d', action='append', default=None,
                        dest='files_dir', required=False, type=str,
                        help="directories with fast5 files to train on")
    parser.add_argument('--ref', '-r', action='store', default=None,
                        dest='ref', required=False, type=str,
                        help="location of refrerence sequence in FASTA")
    parser.add_argument('--iterations', '-i', action='store', dest='iter', default=10,
                        required=False, type=int)
    parser.add_argument('--train_amount', '-a', action='store', dest='amount', default=15,
                        required=False, type=int,
                        help="limit the total length of sequence to use in training.")
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False, default=4,
                        type=int, help="number of jobs to run concurrently")
    parser.add_argument('--cytosine_substitution', '-cs', action='append', default=None,
                        dest='cytosine_sub', required=False, type=str,
                        help="mutate cytosines to this letter in the reference")

    parser.add_argument('--out', '-o', action='store', type=str, required=True, dest='out')

    return parser.parse_args()


def get_hdp_type(requested_type):
        hdp_types = {
            "singleLevelFixed": 0,
            "singleLevelPrior": 1,
            "multisetFixed": 2,
            "multisetPrior": 3,
        }
        assert (requested_type in hdp_types.keys()), "Requested HDP type is invalid, got {}".format(requested_type)
        return hdp_types[requested_type]


def get_initial_hdp_args(args, hdp_type):
    if hdp_type == 0 or hdp_type == 2:
        assert None not in [args.base_gamma, args.leaf_gamma], \
            "ERROR: need to specify concentration parameters for type {}".format(hdp_type)
        if hdp_type == 0:
            return "-B {base} -L {leaf} ".format(base=args.base_gamma, leaf=args.leaf_gamma)
        if hdp_type == 2:
            assert args.middle_gamma is not None, "ERROR: need to specify middle concentration param"
            return "-B {base} -M {middle} -L {leaf} ".format(base=args.base_gamma, middle=args.middle_gamma,
                                                             leaf=args.leaf_gamma)
    else:
        assert None not in [args.base_alpha, args.base_beta, args.leaf_alpha, args.leaf_beta], \
                "ERROR: missing Gamma prior hyper parameters"
        if hdp_type == 1:
            return "-g {Ba} -r {Bb} -i {La} -u {Lb} ".format(Ba=args.base_alpha, Bb=args.base_beta,
                                                             La=args.leaf_alpha, Lb=args.leaf_beta)
        if hdp_type == 3:
            assert None not in [args.middle_alpha, args.middle_beta], "ERROR: need middle hyper parameters"
            return "-g {Ba} -r {Bb} -j {Ma} -y {Mb} -i {La} -u {Lb} ".format(Ba=args.base_alpha, Bb=args.base_beta,
                                                                             Ma=args.middle_alpha, Mb=args.middle_beta,
                                                                             La=args.leaf_alpha, Lb=args.leaf_beta)


# Pipeline Script
args = parse_args()  # parse arguments
working_directory = args.out  # this is the directory we will use for everything
assert os.path.isdir(working_directory), "ERROR: the working directory you specified doesn't exist."
pipeline_log = open(working_directory + "pipeline.log", 'a')
command_line = " ".join(sys.argv[:])
pipeline_log.write("[pipeline] Command Line: {}\n".format(command_line))

# build alignment
signalAlign_directory = "../../signalAlign/"
build_alignment_location = working_directory + "buildAlignment.tsv"
build_alignment_command = "{sA}scripts/makeBuildAlignments.py -o={bA} -t={threshold} " \
                          "".format(sA=signalAlign_directory, C=args.C_alns, mC=args.mC_alns, threshold=args.threshold,
                                    hmC=args.hmC_alns, bA=build_alignment_location)
if args.C_alns is not None: build_alignment_command += "-C={C} ".format(C=args.C_alns)
if args.mC_alns is not None: build_alignment_command += "-mC={mC} ".format(mC=args.mC_alns)
if args.hmC_alns is not None: build_alignment_command += "-hmC={hmC} ".format(hmC=args.hmC_alns)
pipeline_log.write("[pipeline] NOTICE: Making build alignment using files from:\n\t{C}\n\t{mC}\n\t{hmC}\n"
                   "".format(C=args.C_alns, mC=args.mC_alns, hmC=args.hmC_alns))
pipeline_log.write("[pipeline] Command: {}".format(build_alignment_command))
check_call(build_alignment_command.split(), stderr=pipeline_log, stdout=pipeline_log)

# initial HDP
assert (os.path.isfile(build_alignment_location)), "ERROR: Didn't find build alignment"
assert (os.path.exists("./buildHdpUtil")), "ERROR: Didn't find buildHdpUtil"
pipeline_log.write("[pipeline] NOTICE: Making initial HDP of type {}\n".format(args.hdp_type))
template_hdp_location = working_directory + "template." + args.hdp_type + ".nhdp"
complement_hdp_location = working_directory + "complement." + args.hdp_type + ".nhdp"
initial_hdp_build_out = open(working_directory + "build_initial_hdp.out", 'w')
initial_hdp_build_err = open(working_directory + "build_initial_hdp.err", 'w')
if args.verbose is True:
    verbose_flag = "-o "
else:
    verbose_flag = ""
build_initial_hdp_command = "./buildHdpUtil {verbose}-p {hdpType} -v {tHdpLoc} -w {cHdpLoc} -l {buildAln} " \
                            "-n {samples} -I {burnIn} -t {thin} -s {start} -e {end} -k {len} " \
                            "".format(hdpType=get_hdp_type(args.hdp_type), tHdpLoc=template_hdp_location,
                                      cHdpLoc=complement_hdp_location, buildAln=build_alignment_location,
                                      samples=args.gibbs_samples, burnIn=args.burnIn, thin=args.thinning,
                                      start=args.grid_start, end=args.grid_end, len=args.grid_length,
                                      verbose=verbose_flag)
build_initial_hdp_command += get_initial_hdp_args(args=args, hdp_type=get_hdp_type(args.hdp_type))
pipeline_log.write("[pipeline] Command: {}".format(build_initial_hdp_command))
check_call(build_initial_hdp_command.split(), stdout=initial_hdp_build_out, stderr=initial_hdp_build_err)
initial_hdp_build_out.close()
initial_hdp_build_err.close()

# trainModels
assert (os.path.isfile(template_hdp_location) and os.path.isfile(complement_hdp_location)), "ERROR: couldn't find HDPs"
pipeline_log.write("[pipeline] NOTICE: Training HDP models.\n")
template_trained_hdp_location = working_directory + "template_trained." + args.hdp_type + ".nhdp"
complement_trained_hdp_location = working_directory + "complement_trained." + args.hdp_type + ".nhdp"
# make a copy of the HDP files so we can compare before and after training
copyfile(template_hdp_location, template_trained_hdp_location)
copyfile(complement_hdp_location, complement_trained_hdp_location)
train_hdp_out = open(working_directory + "train_hdp.out", 'w')
train_hdp_err = open(working_directory + "train_hdp.err", 'w')
train_models_command = "./trainModels -r={ref} -i={iter} -a={amount} -smt=threeStateHdp -tH={tHdp} " \
                       "-cH={cHdp} -o={wd} -t={threshold} -s={samples} -I={burnIn} -th={thinning} " \
                       "".format(ref=args.ref, iter=args.iter, amount=args.amount, tHdp=template_trained_hdp_location,
                                 cHdp=complement_trained_hdp_location, wd=working_directory, threshold=args.threshold,
                                 samples=args.gibbs_samples, burnIn=args.burnIn, thinning=args.thinning)
assert (len(args.files_dir) >= 1), "ERROR: need to provide at least 1 directory of reads to train on."
for directory in args.files_dir:
    train_models_command += "-d={dir} ".format(dir=directory)
if args.cytosine_sub is not None:
    pipeline_log.write("[pipeline] NOTICE: using cytosine substitutions\n")
    # TODO fill with None?
    assert len(args.cytosine_sub) == len(args.files_dir), "ERROR: need to provide a cytosine substitution for each " \
                                                          "directory.  Just use C if you don't want a change."
    for substitution in args.cytosine_sub:
        train_models_command += "-cs={sub} ".format(sub=substitution)
pipeline_log.write("[pipeline] Command: {}".format(train_models_command))
check_call(train_models_command.split(), stdout=train_hdp_out, stderr=train_hdp_err)

# get HDP distributions
pipeline_log.write("[pipeline] NOTICE: running compareDistributions.\n")
template_trained_distr_dir = working_directory + "template_distrs/"
complement_trained_distr_dir = working_directory + "complement_distrs/"
template_untrained_distr_dir = working_directory + "template_distrs_untrained/"
complement_untrained_distr_dir = working_directory + "complement_distrs_untrained/"
os.makedirs(template_trained_distr_dir)
os.makedirs(complement_trained_distr_dir)
compare_distributions_commands = [
    "./compareDistributions {tHdp} {tDir}".format(tHdp=template_trained_hdp_location,
                                                  tDir=template_trained_distr_dir),
    "./compareDistributions {cHdp} {cDir}".format(cHdp=complement_trained_hdp_location,
                                                  cDir=complement_trained_distr_dir),
    "./compareDistributions {tHdp} {tDir}".format(tHdp=template_hdp_location,
                                                  tDir=template_untrained_distr_dir),
    "./compareDistributions {cHdp} {cDir}".format(cHdp=complement_hdp_location,
                                                  cDir=complement_untrained_distr_dir)
]
for command in compare_distributions_commands:
    pipeline_log.write("[pipeline] Command {}\n".format(command))
procs = [Popen(x.split(), stdout=pipeline_log, stderr=pipeline_log) for x in compare_distributions_commands]
status = [p.wait() for p in procs]
pipeline_log.write("\n[pipeline] DONE.\n")
pipeline_log.close()





