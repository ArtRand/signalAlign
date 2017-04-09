#!/usr/bin/env python
"""Train HMMs for alignment of signal data from the MinION
"""
from __future__ import print_function, division

import sys
import os
import urlparse
import textwrap
import yaml
import h5py

from argparse import ArgumentParser
from random import shuffle
from shutil import copyfile

from multiprocessing import Process, current_process, Manager

from signalalign import parseFofn, DEFAULT_TRAINMODELS_OPTIONS
from signalalign.signalAlignment import SignalAlignment
from signalalign.hiddenMarkovModel import ContinuousPairHmm, HdpSignalHmm
from signalalign.utils import processReferenceFasta
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.utils.bwaWrapper import getBwaIndex


class AbstractSamples(object):
    def __init__(self, source, reference_map):
        self.source = source
        self.reference_map = reference_map

    def _parse(self):
        raise NotImplementedError

    def getFiles(self):
        raise NotImplementedError

    def getKey(self):
        return self.source

    def getReferenceMap(self):
        return self.reference_map


class Fast5Directory(AbstractSamples):
    def __init__(self, source, reference_map):
        AbstractSamples.__init__(self, source, reference_map)
        self.files = self._parse()

    def _parse(self):
        return [self.source + x for x in os.listdir(self.source) if x.endswith(".fast5")]

    def getFiles(self):
        return self.files


class FileOfFilenames(AbstractSamples):
    def __init__(self, source, reference_map):
        AbstractSamples.__init__(self, source, reference_map)
        self.files = self._parse()

    def _parse(self):
        return parseFofn(self.source)

    def getFiles(self):
        return self.files


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
                        dest='substitution_file', help="Substitution positions")
    parser.add_argument("--motif", action="store", dest="motif_key", default=None)
    parser.add_argument('--ambig_char', '-X', action='append', required=False, default=None, type=str, dest='labels',
                        help="Character to substitute at positions, default is 'X'.")
    parser.add_argument('--diagonalExpansion', '-e', action='store', dest='diag_expansion', type=int,
                        required=False, default=None,
                        help="number of diagonals to expand around each anchor default: 50")
    parser.add_argument('--constraintTrim', '-m', action='store', dest='constraint_trim', type=int,
                        required=False, default=None, help='amount to remove from an anchor constraint')
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


def cull_training_files(samples, training_amount, twoD):
    print("trainModels - culling training files.\n", end="", file=sys.stderr)
    training_files = []
    for sample in samples:
        shuffle(sample.getFiles())
        total_amount    = 0
        file_count      = 0
        get_seq_len_fcn = get_2d_length if twoD else get_1d_length
        # loop over files and add them to training list, break when we have enough bases to complete a batch
        # make a list of tuples [(fast5_path, (plus_ref_seq, minus_ref_seq))]
        for f in sample.getFiles():
            training_files.append((f, sample.getReferenceMap()))
            file_count += 1
            total_amount += get_seq_len_fcn(f)
            if total_amount >= training_amount:
                break
        print("Culled {file_count} training files, for {bases} from {sample}."
              .format(file_count=file_count, bases=total_amount, sample=sample.getKey()), end="\n", file=sys.stderr)

    shuffle(training_files)
    return training_files  # [(path_to_fast5, reference_map)...]


def get_expectations(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            alignment = SignalAlignment(**f)
            alignment.run(get_expectations=True)
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def get_model(model_type, model_file):
    assert (model_type in ["threeState", "threeStateHdp"]), "Unsupported StateMachine type"
    # todo clean this up
    if model_type == "threeState":
        assert model_file is not None, "Need to have starting lookup table for {} HMM".format(type)
        model = ContinuousPairHmm(model_type=model_type)
        model.load_model(model_file=model_file)
        return model
    if model_type == "threeStateHdp":
        model = HdpSignalHmm(model_type=model_type)
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


def validateConfig(config):
    # check for inputs
    if config["fast5_dir"] is None and config["fofn"] is None:
        raise RuntimeError("Need to provide a directory of Fast5 files or a file of filenames (fofn)")

    # check for valid paths (if local)
    ref_url = urlparse.urlparse(config["reference_url"])
    if ref_url.scheme == "file":
        if not os.path.exists(ref_url.path):
            raise RuntimeError("Cannot find file: %s" % config["reference_url"])
    return


def generateConfig(config_path):
        if os.path.exists(config_path):
            raise RuntimeError
        config_content = textwrap.dedent("""\
                # SignalAlign model training config file
                output_dir: ../tests/
                samples: [
                {
                    fast5_dir:  ../tests/minion_test_reads/C/,
                    fofn:,
                    positions_file:,
                    motif:,
                    label:,
                }
                ]
                reference:  ../tests/test_sequences/zymo_sequence.fasta
                bwt:
                stateMachineType: threeState
                in_T_Hmm: ../models/testModelR73_acegot_template.model
                in_C_Hmm: ../models/testModelR73_acegot_complement.model
                templateHdp:
                complementHdp:
                iterations: 2
                training_bases: 1000
                job_count: 4
                diagonal_expansion:
                constraint_trim:
                twoD: true

                DEBUG: true
                TEST:

                """)
        fH = open(config_path, "w")
        fH.write(config_content)
        fH.flush()
        fH.close()


def trainModelTransitions(config):
    def process_sample(sample):
        options = dict(**DEFAULT_TRAINMODELS_OPTIONS)
        options.update(sample)
        if options["fast5_dir"] is None and options["fofn"] is None:
            raise RuntimeError("Need to provide path to .fast5 files or file with filenames (fofn)")
        reference_map = processReferenceFasta(fasta=config["reference"],
                                              work_folder=working_folder,
                                              motif_key=options["motif"],
                                              sub_char=options["label"],
                                              positions_file=options["positions_file"])
        if options["fast5_dir"] is not None:
            if options["fofn"] is not None:
                print("WARNING Only using files is directory %s ignoring fofn %s"
                      % (options["files_dir"], options["fofn"]))
            sample = Fast5Directory(options["fast5_dir"], reference_map)
        else:
            sample = FileOfFilenames(options["fofn"], reference_map)
        return sample

    # make directory to put the files we're using
    working_folder = FolderHandler()
    working_folder_path = working_folder.open_folder(config["output_dir"] + "temp_trainModels")
    samples = [process_sample(s) for s in config["samples"]]

    if config["bwt"] is not None:
        print("[trainModels]Using provided BWT")
        bwa_ref_index = config["bwt"]
    else:
        print("signalAlign - indexing reference", file=sys.stderr)
        bwa_ref_index = getBwaIndex(config["reference"], working_folder_path)
        print("signalAlign - indexing reference, done", file=sys.stderr)

    template_model_path   = config["in_T_Hmm"]
    complement_model_path = config["in_C_Hmm"]
    assert os.path.exists(template_model_path) and os.path.exists(complement_model_path), \
        "Missing input models %s and %s" % (template_model_path, complement_model_path)
    template_model   = get_model(config["stateMachineType"], template_model_path)
    complement_model = get_model(config["stateMachineType"], complement_model_path) if config["twoD"] else None

    # get the input HDP, if we're using it
    if config["stateMachineType"] == "threeStateHdp":
        template_hdp   = working_folder.add_file_path("%s" % config["templateHdp"].split("/")[-1])
        copyfile(config["templateHdp"], template_hdp)
        if config["twoD"]:
            complement_hdp = working_folder.add_file_path("%s" % config["complementHdp"].split("/")[-1])
            copyfile(config["complementHdp"], complement_hdp)
        else:
            complement_hdp = None
    else:
        template_hdp   = None
        complement_hdp = None

    # make some paths to files to hold the HMMs
    template_hmm     = working_folder.add_file_path("template_trained.hmm")
    complement_hmm   = working_folder.add_file_path("complement_trained.hmm")
    trained_models   = [template_hmm, complement_hmm]
    untrained_models = [template_model_path, complement_model_path]

    for default_model, trained_model in zip(untrained_models, trained_models):
        assert os.path.exists(default_model), "Didn't find default model {}".format(default_model)
        copyfile(default_model, trained_model)
        assert os.path.exists(trained_model), "Problem copying default model to {}".format(trained_model)

    # start iterating
    i = 0
    while i < config["iterations"]:
        # first cull a set of files to get expectations on
        training_files = cull_training_files(samples=samples,
                                             training_amount=config["training_bases"],
                                             twoD=config["twoD"])
        # setup
        workers    = config["job_count"]
        work_queue = Manager().Queue()
        done_queue = Manager().Queue()
        jobs = []

        # get expectations for all the files in the queue
        # file_ref_tuple should be (fast5, (plus_ref_seq, minus_ref_seq))
        for fast5, ref_map in training_files:
            alignment_args = {
                "reference_map": ref_map,
                "destination": working_folder_path,
                "stateMachineType": config["stateMachineType"],
                "bwa_index": bwa_ref_index,
                "in_templateHmm": template_hmm,
                "in_complementHmm": complement_hmm,
                "in_templateHdp": template_hdp,
                "in_complementHdp": complement_hdp,
                "in_fast5": fast5,
                "threshold": 0.01,
                "diagonal_expansion": config["diagonal_expansion"],
                "constraint_trim": config["constraint_trim"],
                "target_regions": None,
                "degenerate": None,
                "twoD_chemistry": config["twoD"],
            }
            if config["DEBUG"]:
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
        template_expectations_files = [x for x in os.listdir(working_folder_path)
                                       if x.endswith(".template.expectations")]

        complement_expectations_files = [x for x in os.listdir(working_folder_path)
                                         if x.endswith(".complement.expectations")]

        if len(template_expectations_files) > 0:
            add_and_norm_expectations(path=working_folder_path,
                                      files=template_expectations_files,
                                      model=template_model,
                                      hmm_file=template_hmm,
                                      update_transitions=True)

        if config["twoD"] and len(complement_expectations_files) > 0:
            add_and_norm_expectations(path=working_folder_path,
                                      files=complement_expectations_files,
                                      model=complement_model,
                                      hmm_file=complement_hmm,
                                      update_transitions=True)

        # log the running likelihood
        if len(template_model.running_likelihoods) > 0 and \
                (config["twoD"] and len(complement_model.running_likelihoods)) > 0:
            print("{i}| {t_likelihood}\t{c_likelihood}".format(t_likelihood=template_model.running_likelihoods[-1],
                                                               c_likelihood=complement_model.running_likelihoods[-1],
                                                               i=i))
            if config["TEST"] and (len(template_model.running_likelihoods) >= 2) and \
                    (config["twoD"] and len(complement_model.running_likelihoods) >= 2):
                print("TESTING")
                assert (template_model.running_likelihoods[-2] < template_model.running_likelihoods[-1]) and \
                       (complement_model.running_likelihoods[-2] < complement_model.running_likelihoods[-1]), \
                    "Testing: Likelihood error, went up"
        i += 1

    # if we're using HDP, trim the final Hmm (remove assignments)

    print("trainModels - finished training routine", file=sys.stdout)
    print("trainModels - finished training routine", file=sys.stderr)


def main():
    def parse_args():
        parser = ArgumentParser()
        subparsers = parser.add_subparsers(dest="command")

        # parsers for running the full pipeline
        run_parser = subparsers.add_parser("run", help="runs full workflow ")
        run_parser.add_argument('--config', default='trainModels-config.yaml', type=str,
                                help='Path to the (filled in) config file, generated with "generate".')
        subparsers.add_parser("generate", help="generates a config file for your run, do this first")
        return parser.parse_args()
    args   = parse_args()
    if args.command == "generate":
        try:
            config_path = os.path.join(os.getcwd(), "trainModels-config.yaml")
            generateConfig(config_path)
        except RuntimeError:
            print("Using existing config file {}".format(config_path))
            pass
    elif args.command == "run":
        if not os.path.exists(args.config):
            print("{config} not found run generate-config".format(config=args.config))
            exit(1)
        # Parse config
        config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        trainModelTransitions(config)


if __name__ == "__main__":
    sys.exit(main())
