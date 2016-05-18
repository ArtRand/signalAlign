#!/usr/bin/env python

import sys
import unittest
import glob
import os
import shutil
import filecmp
from subprocess import call


class LibTest(unittest.TestCase):
    def test_signalAlign_library(self):
        command = "./signalAlignLibTests"
        result = call(command, shell=True, bufsize=-1, stdout=sys.stdout, stderr=sys.stderr)
        self.assertTrue(result == 0, "signalAlign Library Tests Fail")


class signalAlign_alignment_test(unittest.TestCase):
    def check_alignments(self, true_alignments, reads, reference, test_directory):
        assert len(glob.glob(reads + "*.fast5")) > 0, "Didn't find zymo test MinION reads"
        assert os.path.isfile(reference), "Didn't find zymo reference sequence"

        os.makedirs(test_directory)

        alignment_command = "./runSignalAlign -d={reads} -r={ref} -smt=threeState -o={testDir} " \
                            "".format(reads=reads, ref=reference, testDir=test_directory)
        null_output = open(os.devnull, 'w')
        result = call(alignment_command, shell=True, bufsize=-1, stdout=null_output, stderr=null_output)

        self.assertTrue(result == 0, "error running signalAlign on zymo reads, command was {}"
                                     "".format(alignment_command))

        test_alignments = glob.glob(test_directory + "tempFiles_alignment/*.tsv")

        for alignment in test_alignments:
            alignment_file = alignment.split("/")[-1]
            self.assertTrue(filecmp.cmp(alignment, true_alignments + alignment_file))
        shutil.rmtree(test_directory)

    def test_zymo_reads(self):
        signalAlign_root_directory = "../"
        zymo_true_alignments = signalAlign_root_directory + "tests/test_alignments/zymo_C_test_alignments_sm3/" \
                                                            "tempFiles_alignment/"
        zymo_c_reads = signalAlign_root_directory + "tests/minion_test_reads/C/"
        zymo_reference = signalAlign_root_directory + "tests/test_sequences/zymo_sequence.fasta"
        self.check_alignments(true_alignments=zymo_true_alignments, reads=zymo_c_reads, reference=zymo_reference,
                              test_directory="./test_zymo_reads/")

    def test_ecoli_reads(self):
        signalAlign_root_directory = "../../signalAlign/"
        ecoli_true_alignments = signalAlign_root_directory + "tests/test_alignments/ecoli_test_alignments_sm3/" \
                                                             "tempFiles_alignment/"
        ecoli_reads = signalAlign_root_directory + "tests/minion_test_reads/ecoli/"
        ecoli_reference = signalAlign_root_directory + "tests/test_sequences/E.coli_K12.fasta"
        self.check_alignments(true_alignments=ecoli_true_alignments, reads=ecoli_reads, reference=ecoli_reference,
                              test_directory="./test_ecoli_reads/")


def main():
    testSuite = unittest.TestSuite()
    testSuite.addTest(LibTest('test_signalAlign_library'))
    testSuite.addTest(signalAlign_alignment_test('test_zymo_reads'))
    #testSuite.addTest(signalAlign_alignment_test('test_ecoli_reads'))

    testRunner = unittest.TextTestRunner(verbosity=1)
    testRunner.run(testSuite)

if __name__ == '__main__':
    main()