#!/usr/bin/env python

from __future__ import print_function
import sys
import os
from random import shuffle
from signalAlignLib import get_npRead_2dseq_and_models, exonerated_bwa
from argparse import ArgumentParser
from serviceCourse.file_handlers import FolderHandler


def parse_args():
    parser = ArgumentParser(description=__doc__)

    parser.add_argument('--file_directory', '-d', action='store',
                        dest='files_dir', required=True, type=str, default=None,
                        help="directory with MinION fast5 reads to align")
    #parser.add_argument('--ref', '-r', action='store',
    #                    dest='ref', required=True, type=str,
    #                    help="reference sequence to align to, in FASTA")
    parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for template events, if you don't want the default")
    parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for complement events, if you don't want the default")

    parser.add_argument('--threshold', '-t', action='store', dest='threshold', type=float, required=False,
                        default=0.01, help="posterior match probability threshold, Default: 0.01")

    #parser.add_argument('--constraintTrim', '-m', action='store', dest='trim', type=int,
    #                    required=False, default=14, help='amount to remove from an anchor constraint')

    #parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
    #                    default=4, type=int, help="number of jobs to run concurrently")
    parser.add_argument('--nb_files', '-n', action='store', dest='nb_files', required=False,
                        default=500, type=int, help="maximum number of reads to align")

    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the alignments")

    args = parser.parse_args()
    return args


def estimate_params_with_anchors(fast5, working_folder, bwa_index, forward_reference_path, backward_reference_path, threshold):
    assert (isinstance(working_folder, FolderHandler))

    read_name = fast5.split("/")[-1][:-6]  # get the name without the '.fast5'

    npRead_path = working_folder.add_file_path(read_name + ".npRead")
    npRead_fasta = working_folder.add_file_path(read_name + ".2dSeq.fasta")

    success, def_template_model, def_complement_model = get_npRead_2dseq_and_models(fast5=fast5,
                                                                                    npRead_path=npRead_path,
                                                                                    twod_read_path=npRead_fasta)
    if success is False:
        return False

    # get orientation and cigar from BWA this serves as the guide alignment
    cigar_string, strand = exonerated_bwa(bwa_index=bwa_index, query=npRead_fasta)

    if strand == "-":
        return False

    # input (match) models
    template_lookup_table = "../models/testModel_template.model"
    complement_lookup_table = "../models/testModel_complement.model" if def_complement_model else \
        "../models/testModel_complement_pop1.model"

    binary = "./estimateNanoporeParams"

    command = "echo {cigar} | {bin} -T {tLuT} -C {cLuT} -q {npRead} -f {fRef} -b {bRef} -D {threshold}" \
              "".format(cigar=cigar_string, bin=binary, tLuT=template_lookup_table, cLuT=complement_lookup_table,
                        npRead=npRead_path, fRef=forward_reference_path, bRef=backward_reference_path,
                        threshold=threshold)

    print("running command {command}".format(command=command), file=sys.stderr)

    os.system(command)


def estimate_params(fast5, working_folder, threshold):
    assert (isinstance(working_folder, FolderHandler))

    read_name = fast5.split("/")[-1][:-6]  # get the name without the '.fast5'

    npRead_path = working_folder.add_file_path(read_name + ".npRead")
    npRead_fasta = working_folder.add_file_path(read_name + ".2dSeq.fasta")

    success, def_template_model, def_complement_model = get_npRead_2dseq_and_models(fast5=fast5,
                                                                                    npRead_path=npRead_path,
                                                                                    twod_read_path=npRead_fasta)
    if success is False:
        return False

    # input (match) models
    template_lookup_table = "../models/testModelR9_template.model"
    complement_lookup_table = "../models/testModelR9_complement_pop2.model"

    binary = "./estimateNanoporeParams"

    command = "{bin} -T {tLuT} -C {cLuT} -q {npRead} -D {threshold}" \
              "".format(bin=binary, tLuT=template_lookup_table, cLuT=complement_lookup_table, npRead=npRead_path,
                        threshold=threshold)

    print("running command {command}".format(command=command), file=sys.stderr)

    os.system(command)


def main(args):
    args = parse_args()

    # make directory to put temporary files
    temp_folder = FolderHandler()
    temp_dir_path = temp_folder.open_folder(args.out + "npParamEstimation")

    fast5s = [x for x in os.listdir(args.files_dir) if x.endswith(".fast5")]

    if len(fast5s) > args.nb_files:
        shuffle(fast5s)
        fast5s = fast5s[:args.nb_files]
    print("kmer_index\tevent_mean\tevent_stdv\tstart\tp_model\tscale\tshift\tdrift\tstrand\n", file=sys.stdout)
    for fast5 in fast5s:
        #estimate_params(fast5=args.files_dir + fast5, working_folder=temp_folder, bwa_index=bwa_ref_index,
        #                forward_reference_path=plus_strand_sequence, backward_reference_path=minus_strand_sequence,
        #                threshold=args.threshold)
        estimate_params(fast5=args.files_dir + fast5, working_folder=temp_folder, threshold=args.threshold)

    temp_folder.remove_folder()
    return

if __name__ == "__main__":
    sys.exit(main(sys.argv))
