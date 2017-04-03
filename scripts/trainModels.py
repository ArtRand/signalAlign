#!/usr/bin/env python
"""Train HMMs for alignment of signal data from the MinION
"""
from __future__ import print_function, division

import sys
import os
import h5py

from argparse import ArgumentParser
from random import shuffle
from shutil import copyfile

from multiprocessing import Process, current_process, Manager

from signalalign.SignalAlignment import SignalAlignment
from signalalign.hiddenMarkovModel import ContinuousPairHmm, HdpSignalHmm
from signalalign.utils import processReferenceFasta
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.utils.bwaWrapper import getBwaIndex


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--file_directory', '-d', action='append', default=None,
                        dest='files_dir', required=True, type=str,
                        help="directories with fast5 files to train on. example: ../reads/")
    parser.add_argument('--ref', '-r', action='store', default=None, dest='ref', required=True, type=str,
                        help="location of refrerence sequence in FASTA, example: ../ref.fasta")
    parser.add_argument('--output_location', '-o', action='store', dest='out', default=None,
                        required=True, type=str,
                        help="directory to put the trained model, and use for working directory. example: ./scratch/")
    # optional arguments
    parser.add_argument("--2d", action='store_true', dest="twoD", default=False, help="flag, reads are 2D chemistry.")
    parser.add_argument("--bwt", action='store', dest="bwt", default=None,
                        help="path to BWT files. example: ../ref.fasta")
    parser.add_argument('--stateMachineType', '-smt', action='store', dest='stateMachineType', type=str,
                        default="threeState", required=False,
                        help="StateMachine options: threeState, threeStateHdp")
    parser.add_argument("--file_of_files", "-fofn", action="append", required=False, default=None, dest="fofn",
                        type=str, help="text file with absolute paths of files to use")
    parser.add_argument('--iterations', '-i', action='store', dest='iter', default=10,
                        required=False, type=int, help='number of iterations to perform')
    parser.add_argument('--train_amount', '-a', action='store', dest='amount', default=15000,
                        required=False, type=int,
                        help="limit the total length of sequence to use in training (batch size).")
    parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm',
                        required=True, type=str, help="template model to bootstrap from, find a starting model in the "
                        "models directory")
    parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm',
                        required=True, type=str, help="complement model to bootstrap from, find a starting model in the "
                        "models directory")
    parser.add_argument('--templateHDP', '-tH', action='store', dest='templateHDP', default=None,
                        help="path to template HDP model to use")
    parser.add_argument('--complementHDP', '-cH', action='store', dest='complementHDP', default=None,
                        help="path to complement HDP model to use")
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False, default=4,
                        type=int, help="number of jobs to run concurrently")
    parser.add_argument('--test', action='store_true', default=False, dest='test', help="Used for CI testing")
    parser.add_argument('--ambiguity_positions', '-p', action='store', required=False, default=None,
                        dest='substitution_file', help="Ambiguity positions")
    parser.add_argument("--motif", action="store", dest="motif_key", default=None)
    parser.add_argument('--ambig_char', '-X', action='append', required=False, default=None, type=str, dest='ambig_char',
                        help="Character to substitute at positions, default is 'X'.")
    parser.add_argument('--debug', action='store_true', dest="DEBUG", default=False)

    args = parser.parse_args()
    return args


def get_2d_length(fast5):
    read = h5py.File(fast5, 'r')
    read_length = 0
    twoD_read_sequence_address = "/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"
    if not (twoD_read_sequence_address in read):
        print("This read didn't have a 2D read", fast5, end='\n', file=sys.stderr)
        read.close()
        return 0
    else:
        read_length = len(read[twoD_read_sequence_address][()].split()[2])
        read.close()
        return read_length


def get_1d_length(fast5):
    read = h5py.File(fast5, "r")
    read_length = 0
    template_fastq_address = "/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"
    if not (template_fastq_address in read):
        print("Read %s has not been basecalled" % fast5)
        read.close()
        return 0
    else:
        read_length = len(read[template_fastq_address][()].split()[2])
        print("read %s has %s bases" % (fast5, read_length))
        read.close()
        return read_length


def cull_training_files(directories, fofns, training_amount, reference_maps, twoD):
    def get_file_list():
        source_list_tups = []
        if fofns is not None:
            for fofn in fofns:
                source_list_tups.append((fofn, parse_fofn(fofn)))
            return source_list_tups
        else:
            assert directories is not None
            for d in directories:
                fast5s = [d + x for x in os.listdir(d) if x.endswith(".fast5")]
                source_list_tups.append((d, fast5s))
            return source_list_tups

    print("trainModels - culling training files.\n", end="", file=sys.stderr)

    training_files = []
    add_to_training_files = training_files.append

    sources_and_files = get_file_list()

    assert len(sources_and_files) == len(reference_maps), "Need to have equal number of references as training " \
                                                          "file directories."
    # loop over the directories and collect reads for training
    for j, tup in enumerate(sources_and_files):
        assert len(tup) == 2
        source = tup[0]  # the directory or the fofn
        fast5s = tup[1]  # the list of files (paths)
        shuffle(fast5s)
        total_amount = 0
        n = 0
        get_seq_len_fcn = get_2d_length if twoD else get_1d_length
        # loop over files and add them to training list, break when we have enough bases to complete a batch
        # make a list of tuples [(fast5_path, (plus_ref_seq, minus_ref_seq))]
        for i in xrange(len(fast5s)):
            add_to_training_files((fast5s[i], reference_maps[j]))
            n += 1
            total_amount += get_seq_len_fcn(fast5s[i])
            if total_amount >= training_amount:
                break
        print("Culled {nb_files} training files, for {bases} from {dir}.".format(nb_files=n, bases=total_amount,
                                                                                 dir=source),
              end="\n", file=sys.stderr)

    shuffle(training_files)
    return training_files


def get_expectations(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            alignment = SignalAlignment(**f)
            alignment.run(get_expectations=True)
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def get_model(type, threshold, model_file):
    assert (type in ["threeState", "threeStateHdp"]), "Unsupported StateMachine type"
    # todo clean this up
    if type == "threeState":
        assert model_file is not None, "Need to have starting lookup table for {} HMM".format(type)
        model = ContinuousPairHmm(model_type=type)
        model.load_model(model_file=model_file)
        return model
    if type == "threeStateHdp":
        model = HdpSignalHmm(model_type=type, threshold=threshold)
        model.load_model(model_file=model_file)
        return model


def add_and_norm_expectations(path, files, model, hmm_file, update_transitions=False, update_emissions=False):
    if update_emissions is False and update_transitions is False:
        print("[trainModels] NOTICE: Training transitions by default\n", file=sys.stderr)
        update_transitions = True

    model.likelihood = 0
    files_added_successfully = 0
    files_with_problems = 0
    for f in files:
        try:
            success = model.add_expectations_file(path + f)
            os.remove(path + f)
            if success:
                files_added_successfully += 1
            else:
                files_with_problems += 1
        except Exception as e:
            print("Problem adding expectations file {file} got error {e}".format(file=path + f, e=e),
                  file=sys.stderr)
            os.remove(path + f)
            files_with_problems += 1
    model.normalize(update_transitions=update_transitions, update_emissions=update_emissions)
    model.write(hmm_file)
    model.running_likelihoods.append(model.likelihood)
    if type(model) is HdpSignalHmm:
        model.reset_assignments()
    print("[trainModels] NOTICE: Added {success} expectations files successfully, {problem} files had problems\n"
          "".format(success=files_added_successfully, problem=files_with_problems), file=sys.stderr)


def build_hdp(template_hdp_path, complement_hdp_path, template_assignments, complement_assignments, samples,
              burn_in, thinning, verbose=False):
    assert (template_assignments is not None) and (complement_assignments is not None), \
        "trainModels - ERROR: missing assignments"

    if verbose is True:
        verbose_flag = "--verbose "
    else:
        verbose_flag = ""

    command = "./buildHdpUtil {verbose}-v {tHdpP} -w {cHdpP} -E {tExpectations} -W {cExpectations} " \
              "-n {samples} -I {burnIn} -t {thinning}".format(tHdpP=template_hdp_path,
                                                              cHdpP=complement_hdp_path,
                                                              tExpectations=template_assignments,
                                                              cExpectations=complement_assignments,
                                                              samples=samples, burnIn=burn_in,
                                                              thinning=thinning,
                                                              verbose=verbose_flag)
    print("[trainModels] Running command:{}".format(command), file=sys.stderr)
    os.system(command)  # todo try checkoutput
    print("trainModels - built HDP.", file=sys.stderr)
    return


def main(args):
    # parse command line arguments
    args = parse_args()

    command_line = " ".join(sys.argv[:])
    print("Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)

    start_message = """\n
    # Starting Baum-Welch training.
    # Directories with training files: {files_dir}
    # Each batch has {amount} bases, performing {iter} iterations.
    # Using reference sequence: {ref}
    # Writing trained models to: {outLoc}
    # Performing {iterations} iterations.
    # Using model: {model}
    # Using HDPs: {thdp} / {chdp}
    # Training emissions: {emissions}
    #        transitions: {transitions}
    \n
    """.format(files_dir=args.files_dir, amount=args.amount, ref=args.ref, outLoc=args.out, iter=args.iter,
               iterations=args.iter, model=args.stateMachineType, thdp=args.templateHDP, chdp=args.complementHDP,
               emissions=args.emissions, transitions=args.transitions)

    if args.files_dir is None and args.fofn is None:
        print("Need to provide directory with .fast5 files of file of file names", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(args.ref):
        print("Did not find valid reference file", file=sys.stderr)
        sys.exit(1)

    print(start_message, file=sys.stdout)

    # make directory to put the files we're using files
    working_folder = FolderHandler()
    working_directory_path = working_folder.open_folder(args.out + "tempFiles_expectations")

    # if we are performing supervised training with multiple kinds of substitutions, then we
    # need to make a reference sequence for each one
    if args.ambig_char is not None:
        reference_maps = []
        for sub_char in args.ambig_char:
            reference_maps.append(processReferenceFasta(args.ref, args.motif_key, working_folder, sub_char))
    else:
        reference_map = processReferenceFasta(fasta=args.ref, work_folder=working_folder)
        reference_maps = [reference_map]

    # index the reference for bwa
    if args.bwt is not None:
        print("[trainModels]Using provided BWT")
        bwa_ref_index = args.bwt
    else:
        print("signalAlign - indexing reference", file=sys.stderr)
        bwa_ref_index = getBwaIndex(args.ref, working_directory_path)
        print("signalAlign - indexing reference, done", file=sys.stderr)

    template_model_path   = args.in_T_Hmm
    complement_model_path = args.in_C_Hmm
    assert os.path.exists(template_model_path) and os.path.exists(complement_model_path), \
        "Missing default lookup tables"
    # make the model objects, for the threeState model, we're going to parse the lookup table or the premade
    # model, for the HDP model, we just load the transitions
    template_model = get_model(type=args.stateMachineType, threshold=args.threshold, model_file=template_model_path)
    complement_model = get_model(type=args.stateMachineType, threshold=args.threshold, model_file=complement_model_path)

    # get the input HDP, if we're using it
    if args.stateMachineType == "threeStateHdp":
        assert (args.templateHDP is not None) and (args.complementHDP is not None), \
            "Need to provide serialized HDP files for this stateMachineType"
        assert (os.path.isfile(args.templateHDP)) and (os.path.isfile(args.complementHDP)),\
            "Could not find the HDP files"
        template_hdp = working_folder.add_file_path("{}".format(args.templateHDP.split("/")[-1]))
        complement_hdp = working_folder.add_file_path("{}".format(args.complementHDP.split("/")[-1]))
        copyfile(args.templateHDP, template_hdp)
        copyfile(args.complementHDP, complement_hdp)
    else:
        template_hdp = None
        complement_hdp = None

    # make some paths to files to hold the HMMs
    template_hmm = working_folder.add_file_path("template_trained.hmm")
    complement_hmm = working_folder.add_file_path("complement_trained.hmm")

    trained_models = [template_hmm, complement_hmm]

    untrained_models = [template_model_path, complement_model_path]

    for default_model, trained_model in zip(untrained_models, trained_models):
        assert os.path.exists(default_model), "Didn't find default model {}".format(default_model)
        copyfile(default_model, trained_model)
        assert os.path.exists(trained_model), "Problem copying default model to {}".format(trained_model)

    print("Starting {iterations} iterations.\n\n\t    Running likelihoods\ni\tTemplate\tComplement".format(
        iterations=args.iter), file=sys.stdout)

    # start iterating
    i = 0
    while i < args.iter:
        # first cull a set of files to get expectations on
        training_files = cull_training_files(directories=args.files_dir, fofns=args.fofn, training_amount=args.amount,
                                             reference_maps=reference_maps, twoD=args.twoD)
        # setup
        workers = args.nb_jobs
        work_queue = Manager().Queue()
        done_queue = Manager().Queue()
        jobs = []

        # get expectations for all the files in the queue
        # file_ref_tuple should be (fast5, (plus_ref_seq, minus_ref_seq))
        for file_ref_tuple in training_files:
            alignment_args = {
                "reference_map": file_ref_tuple[1],
                "path_to_EC_refs": None,
                "destination": working_directory_path,
                "stateMachineType": args.stateMachineType,
                "bwa_index": bwa_ref_index,
                "in_templateHmm": template_hmm,
                "in_complementHmm": complement_hmm,
                "in_templateHdp": template_hdp,
                "in_complementHdp": complement_hdp,
                "in_fast5": file_ref_tuple[0],  # fast5
                "threshold": args.threshold,
                "diagonal_expansion": args.diag_expansion,
                "constraint_trim": args.constraint_trim,
                "target_regions": None,
                "degenerate": None,
                "twoD_chemistry": args.twoD,
            }
            if args.DEBUG:
                alignment = SignalAlignment(**alignment_args)
                alignment.run(get_expectations=True)
            else:
                work_queue.put(alignment_args)

        for w in xrange(workers):
            p = Process(target=get_expectations, args=(work_queue, done_queue))
            p.start()
            jobs.append(p)
            work_queue.put('STOP')

        for p in jobs:
            p.join()

        done_queue.put('STOP')

        # load then normalize the expectations
        template_expectations_files = [x for x in os.listdir(working_directory_path)
                                       if x.endswith(".template.expectations")]

        complement_expectations_files = [x for x in os.listdir(working_directory_path)
                                         if x.endswith(".complement.expectations")]

        if len(template_expectations_files) > 0:
            add_and_norm_expectations(path=working_directory_path,
                                      files=template_expectations_files,
                                      model=template_model,
                                      hmm_file=template_hmm,
                                      update_emissions=args.emissions,
                                      update_transitions=args.transitions)

        if len(complement_expectations_files) > 0:
            add_and_norm_expectations(path=working_directory_path,
                                      files=complement_expectations_files,
                                      model=complement_model,
                                      hmm_file=complement_hmm,
                                      update_emissions=args.emissions,
                                      update_transitions=args.transitions)

        # Build HDP from last round of assignments
        if args.stateMachineType == "threeStateHdp" and args.emissions is True:
            assert isinstance(template_model, HdpSignalHmm) and isinstance(complement_model, HdpSignalHmm)
            if min(template_model.assignments_record[-1],
                   complement_model.assignments_record[-1]) < args.min_assignments:
                print("[trainModels] not enough assignments at iteration {}, continuing...".format(i),
                      file=sys.stderr)
                i -= 1
                pass
            else:
                total_assignments = max(template_model.assignments_record[-1], complement_model.assignments_record[-1])

                build_hdp(template_hdp_path=template_hdp, complement_hdp_path=complement_hdp,
                          template_assignments=template_hmm, complement_assignments=complement_hmm,
                          samples=args.gibbs_samples, thinning=args.thinning, burn_in=30 * total_assignments,
                          verbose=args.verbose)

        # log the running likelihood
        if len(template_model.running_likelihoods) > 0 and len(complement_model.running_likelihoods) > 0:
            print("{i}| {t_likelihood}\t{c_likelihood}".format(t_likelihood=template_model.running_likelihoods[-1],
                                                               c_likelihood=complement_model.running_likelihoods[-1],
                                                               i=i))
            if args.test and (len(template_model.running_likelihoods) >= 2) and \
                    (len(complement_model.running_likelihoods) >= 2):
                assert (template_model.running_likelihoods[-2] < template_model.running_likelihoods[-1]) and \
                       (complement_model.running_likelihoods[-2] < complement_model.running_likelihoods[-1]), \
                    "Testing: Likelihood error, went up"
        i += 1

    # if we're using HDP, trim the final Hmm (remove assignments)

    print("trainModels - finished training routine", file=sys.stdout)
    print("trainModels - finished training routine", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
