#!/usr/bin/env python

import sys
import unittest
import glob
import os
import shutil
import pandas as pd
import numpy as np
from subprocess import call

from margin.utils import getFastaDictionary

SIGNALALIGN_ROOT = "../"
ZYMO_C_READS = SIGNALALIGN_ROOT + "tests/minion_test_reads/C/"
ZYMO_REFERENCE = SIGNALALIGN_ROOT + "tests/test_sequences/zymo_sequence.fasta"


def parse_alignment_full(alignment_file):
    data = pd.read_table(alignment_file, usecols=(1, 2, 4, 5, 9, 12, 13),
                         dtype={'ref_pos': np.int64,
                                'ref_kmer': np.str,
                                'strand': np.str,
                                'event_index': np.int64,
                                'kmer': np.str,
                                'posterior_prob': np.float64,
                                'event_mean': np.float64},
                         header=None,
                         names=['ref_pos', 'ref_kmer', 'strand', 'event_index', 'kmer', 'posterior_prob', 'event_mean'])
    return data


class LibTest(unittest.TestCase):
    def test_signalAlign_library(self):
        command = "./signalAlignLibTests"
        result = call(command, shell=True, bufsize=-1, stdout=sys.stdout, stderr=sys.stderr)
        self.assertTrue(result == 0, "signalAlign Library Tests Fail")


class signalAlignLibTests(unittest.TestCase):
    def setUp(self):
        self.work_dir = "./signalAlign_pylibTest/"
        os.makedirs(self.work_dir)

    def tearDown(self):
        shutil.rmtree(self.work_dir)


class SignalAlignAlignmentTest(unittest.TestCase):
    def setUp(self):
        os.makedirs("./signalAlign_unittest/")

    def tearDown(self):
        shutil.rmtree("./signalAlign_unittest/")

    def check_alignments(self, true_alignments, reads, reference, kmer_length, contig_name, extra_args=None):

        def get_kmer(start):
            return referece_sequence[start:start + kmer_length]

        assert len(glob.glob(reads + "*.fast5")) > 0, "Didn't find zymo test MinION reads"
        assert os.path.isfile(reference), "Didn't find zymo reference sequence"

        alignment_command = "./runSignalAlign -d={reads} -r={ref} -smt=threeState -o={testDir} " \
                            "".format(reads=reads, ref=reference, testDir="./signalAlign_unittest/")
        if extra_args is not None:
            alignment_command += extra_args

        null_output = open(os.devnull, 'w')
        result = call(alignment_command, shell=True, bufsize=-1, stdout=null_output, stderr=null_output)

        self.assertTrue(result == 0, "error running signalAlign alignments command was {}"
                                     "".format(alignment_command))

        test_alignments = glob.glob("./signalAlign_unittest/tempFiles_alignment/*.tsv")

        referece_sequence = getFastaDictionary(reference)[contig_name]

        self.assertTrue(len(test_alignments) == len(glob.glob(true_alignments + "*.tsv")),
                        "Didn't make all alignments got {got} should be {should}".format(got=len(test_alignments),
                                                                                         should=len(glob.glob(true_alignments + "*.tsv"))))

        for alignment in test_alignments:
            alignment_file = alignment.split("/")[-1]
            expected = parse_alignment_full(true_alignments + alignment_file)
            obs = parse_alignment_full(alignment)
            self.assertTrue(len(obs) == len(expected))
            for row in obs.itertuples():
                ref_pos = row[1]
                obs_kmer = row[2]
                strand = row[3]
                exp_kmer = get_kmer(ref_pos)
                self.assertTrue(obs_kmer == exp_kmer, msg="kmer at index {idx} on strand {strand} is {obs} should be "
                                                          "{exp}, file {f}".format(idx=ref_pos,
                                                                                   strand=strand,
                                                                                   obs=obs_kmer,
                                                                                   exp=exp_kmer,
                                                                                   f=alignment))

    def test_zymo_reads(self):
        zymo_true_alignments = SIGNALALIGN_ROOT + "tests/test_alignments/zymo_C_test_alignments_sm3/" \
                                                  "tempFiles_alignment/"
        self.check_alignments(true_alignments=zymo_true_alignments,
                              reads=ZYMO_C_READS,
                              reference=ZYMO_REFERENCE,
                              kmer_length=6,
                              contig_name="ZYMO",
                              extra_args="--2d ")

    def test_pUC_r9_reads_5mer(self):
        pUC_true_alignments = SIGNALALIGN_ROOT + "tests/test_alignments/pUC_5mer_tempFiles_alignment/"
        pUC_reads = SIGNALALIGN_ROOT + "tests/minion_test_reads/pUC/"
        pUC_reference = SIGNALALIGN_ROOT + "tests/test_sequences/pUC19_SspI.fa"
        self.check_alignments(true_alignments=pUC_true_alignments,
                              reads=pUC_reads,
                              reference=pUC_reference,
                              kmer_length=5,
                              contig_name="pUC19",
                              extra_args="-T=../models/testModelR9_5mer_acegot_template.model "
                                         "-C=../models/testModelR9_5mer_acegot_complement.model "
                                         "--2d ")

    def test_pUC_r9_reads_6mer(self):
        pUC_true_alignments = SIGNALALIGN_ROOT + "tests/test_alignments/pUC_6mer_tempFiles_alignment/"
        pUC_reads = SIGNALALIGN_ROOT + "tests/minion_test_reads/pUC/"
        pUC_reference = SIGNALALIGN_ROOT + "tests/test_sequences/pUC19_SspI.fa"
        self.check_alignments(true_alignments=pUC_true_alignments,
                              reads=pUC_reads,
                              reference=pUC_reference,
                              kmer_length=6,
                              contig_name="pUC19",
                              extra_args="--2d ")

    def test_Ecoli1D_reads_5mer(self):
        pUC_true_alignments = SIGNALALIGN_ROOT + "tests/test_alignments/ecoli1D_test_alignments_sm3/"
        pUC_reads = SIGNALALIGN_ROOT + "tests/minion_test_reads/1D/"
        pUC_reference = SIGNALALIGN_ROOT + "tests/test_sequences/E.coli_K12.fasta"
        self.check_alignments(true_alignments=pUC_true_alignments,
                              reads=pUC_reads,
                              reference=pUC_reference,
                              kmer_length=5,
                              contig_name="gi_ecoli",
                              extra_args="-T=../models/testModelR9p4_5mer_acegt_template.model ")


class signalAlign_EM_test(unittest.TestCase):
    def setUp(self):
        os.makedirs("./signalAlign_unittest/")

    def tearDown(self):
        shutil.rmtree("./signalAlign_unittest/")

    def test_EM(self):
        em_command = "./trainModels run --config=../tests/trainModels-config.yaml"
        null_output = open(os.devnull, 'w')
        result = call(em_command, shell=True, bufsize=-1, stdout=null_output, stderr=null_output)

        self.assertTrue(result == 0, "error running signalAlign alignments command was {}"
                                     "".format(em_command))


def main():
    testSuite = unittest.TestSuite()
    testSuite.addTest(LibTest('test_signalAlign_library'))
    testSuite.addTest(SignalAlignAlignmentTest('test_zymo_reads'))
    testSuite.addTest(SignalAlignAlignmentTest('test_pUC_r9_reads_5mer'))
    testSuite.addTest(SignalAlignAlignmentTest('test_pUC_r9_reads_6mer'))
    testSuite.addTest(SignalAlignAlignmentTest('test_Ecoli1D_reads_5mer'))
    testSuite.addTest(signalAlign_EM_test('test_EM'))

    testRunner = unittest.TextTestRunner(verbosity=2)
    testRunner.run(testSuite)


if __name__ == '__main__':
    main()
