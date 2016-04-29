#!/usr/bin/env python
"""A collection of functions for visualizing kmer distributions
"""
import os
import sys
sys.path.append("../")
from signalAlignLib import TemplateModel, ComplementModel
import string
import numpy as np
from sklearn.neighbors import KernelDensity
from scipy.stats import gaussian_kde
from scipy.stats import norm


def nucleotideToIndex(base):
    if base == 'A':
        return 0
    if base == 'C':
        return 1
    if base == 'E':
        return 2
    if base == 'G':
        return 3
    if base == 'O':
        return 4
    if base == 'T':
        return 5
    if base == 'N':
        return 6


def getKmerIndex(kmer):
    """This is the algorithm for finding the rank (index) of a kmer)
    """
    alphabet = "ACEGOT"
    axisLength = len(alphabet)**len(kmer)
    l = axisLength/len(alphabet)
    i = 0
    index = 0
    while l > 1:
        index += l*nucleotideToIndex(kmer[i])
        i += 1
        l = l/len(alphabet)
    index += nucleotideToIndex(kmer[-1])
    return int(index)


class KmerDistribution(object):
    def __init__(self, data_directory):
        self.data_directory = data_directory


class KmerAlignmentHistogram(KmerDistribution):
    def __init__(self, data_directory, kmer):
        super(KmerAlignmentHistogram, self).__init__(data_directory=data_directory)
        self.kmer = kmer
        self.data_file = data_directory + "{}_hist.txt".format(self.kmer)

        assert (os.path.isfile(self.data_file)), "No data for kmer {kmer}".format(kmer=self.kmer)

        self.histogram = None
        self.x_vals = None
        self.kde_pdf = None

        self.parse_histogram()

    def parse_histogram(self):
        self.histogram = np.loadtxt(self.data_file, dtype=np.float64)
        self.histogram = [x for x in self.histogram if 0 < x < 100]

    def parse_xvals(self, x_vals_file):
        self.x_vals = np.loadtxt(x_vals_file, dtype=np.float64)

    def make_kde(self, x_vals):
        kernel = gaussian_kde(self.histogram)
        KDE_PDF = kernel.evaluate(x_vals)
        return KDE_PDF


class KmerHdpDistribution(KmerDistribution):
    def __init__(self, data_directory, kmer):
        super(KmerHdpDistribution, self).__init__(data_directory=data_directory)
        self.kmer = kmer
        self.data_file = data_directory + "{}_distr.txt".format(self.kmer)

        assert (os.path.isfile(self.data_file)), "Didn't find distribution file at {}".format(self.data_file)

        self.density = None
        self.parse_density_file()

    def parse_density_file(self):
        self.density = np.loadtxt(self.data_file, dtype=np.float64)
    

def get_kmer_densities(path, kmer):
    mC_trans = string.maketrans("C", "E")
    hmC_trans = string.maketrans("C", "O")
    c_density = KmerHdpDistribution(path, kmer)
    mc_density = KmerHdpDistribution(path, kmer.translate(mC_trans))
    hmc_density = KmerHdpDistribution(path, kmer.translate(hmC_trans))
    return c_density, mc_density, hmc_density


def plot_ont_distribution(kmer, fast5, x_vals):
    def get_model_table(model):
        return model.get_model_dict()
    template_model = get_model_table(TemplateModel(fast5File=fast5))
    complement_model = get_model_table(ComplementModel(fast5File=fast5))

    return norm.pdf(x_vals, template_model[kmer][0], template_model[kmer][1]), \
           norm.pdf(x_vals, complement_model[kmer][0], complement_model[kmer][1])


def plot_hmm_distribution(kmer, hmm, x_vals):
    def get_model_from_hmm(model_file):
        fH = open(hmm, 'r')

        # line 0, type, stateNumber, etc
        line = map(float, fH.readline().split())
        assert len(line) == 4, "Bad header line"

        # line 1, transitions
        line = map(float, fH.readline().split())
        assert len(line) == 10, "Bad transitions line"

        # line 2, model
        line = map(float, fH.readline().split())
        assert len(line) == 6**6 * 2  # number of kmers * normal distribution parameters
        return line

    model = get_model_from_hmm(hmm)
    kmer_index = getKmerIndex(kmer)
    table_index = kmer_index * 2
    print model[table_index], model[table_index + 1]
    return norm.pdf(x_vals, model[table_index], model[table_index + 1])

