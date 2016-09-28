#!/usr/bin/env python

import sys
sys.path.append("../")
import unittest
import glob
import os
import shutil
import pandas as pd
import numpy as np
from subprocess import call
from alignmentAnalysisLib import get_first_sequence

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


class SignalAlignAlignmentTest(unittest.TestCase):
    def setUp(self):
        os.makedirs("./signalAlign_unittest/")

    def tearDown(self):
        shutil.rmtree("./signalAlign_unittest/")

    def check_alignments(self, true_alignments, reads, reference, kmer_length):

        def get_kmer(start):
            return referece_sequence[start:start + kmer_length]

        assert len(glob.glob(reads + "*.fast5")) > 0, "Didn't find zymo test MinION reads"
        assert os.path.isfile(reference), "Didn't find zymo reference sequence"

        alignment_command = "./runSignalAlign -d={reads} -r={ref} -smt=threeState -o={testDir} " \
                            "".format(reads=reads, ref=reference, testDir="./signalAlign_unittest/")
        null_output = open(os.devnull, 'w')
        result = call(alignment_command, shell=True, bufsize=-1, stdout=null_output, stderr=null_output)

        self.assertTrue(result == 0, "error running signalAlign alignments command was {}"
                                     "".format(alignment_command))

        test_alignments = glob.glob("./signalAlign_unittest/tempFiles_alignment/*.tsv")

        referece_sequence = get_first_sequence(reference)

        self.assertTrue(len(test_alignments) == len(glob.glob(true_alignments + "*.tsv")),
                        "Didn't make all alignments")

        for alignment in test_alignments:
            alignment_file = alignment.split("/")[-1]
            expected = parse_alignment_full(true_alignments + alignment_file)
            obs = parse_alignment_full(alignment)
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
            self.assertTrue(expected.equals(obs), msg="{} is not the same".format(alignment_file))

    def test_zymo_reads(self):
        zymo_true_alignments = SIGNALALIGN_ROOT + "tests/test_alignments/zymo_C_test_alignments_sm3/" \
                                                            "tempFiles_alignment/"
        self.check_alignments(true_alignments=zymo_true_alignments, reads=ZYMO_C_READS,
                              reference=ZYMO_REFERENCE, kmer_length=6)

    def test_ecoli_reads(self):
        ecoli_true_alignments = SIGNALALIGN_ROOT + "tests/test_alignments/ecoli_test_alignments_sm3/" \
                                                             "tempFiles_alignment/"
        ecoli_reads = SIGNALALIGN_ROOT + "tests/minion_test_reads/ecoli/"
        ecoli_reference = SIGNALALIGN_ROOT + "tests/test_sequences/E.coli_K12.fasta"
        self.check_alignments(true_alignments=ecoli_true_alignments, reads=ecoli_reads,
                              reference=ecoli_reference, kmer_length=6)


class signalAlign_EM_test(unittest.TestCase):
    def setUp(self):
        os.makedirs("./signalAlign_unittest/")

    def tearDown(self):
        shutil.rmtree("./signalAlign_unittest/")

    def test_EM(self):
        em_command = "./trainModels -d={reads} -r={ref} -o={testDir} -i=3 -a=10000 --transitions -smt=threeState " \
                     "--test".format(reads=ZYMO_C_READS, ref=ZYMO_REFERENCE, testDir="./signalAlign_unittest/")
        null_output = open(os.devnull, 'w')
        result = call(em_command, shell=True, bufsize=-1, stdout=null_output, stderr=null_output)

        self.assertTrue(result == 0, "error running signalAlign alignments command was {}"
                                     "".format(em_command))


def main():
    testSuite = unittest.TestSuite()
    testSuite.addTest(LibTest('test_signalAlign_library'))
    testSuite.addTest(SignalAlignAlignmentTest('test_zymo_reads'))
    #testSuite.addTest(signalAlign_alignment_test('test_ecoli_reads'))
    testSuite.addTest(signalAlign_EM_test('test_EM'))

    testRunner = unittest.TextTestRunner(verbosity=1)
    testRunner.run(testSuite)

if __name__ == '__main__':
    main()
