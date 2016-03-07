#!/usr/bin/env python
"""A collection of functions for visualizing kmer distributions
"""
import os
import sys
sys.path.append("../")
from signalAlignLib import TemplateModel, ComplementModel
import string
import numpy as np
from scipy.stats.kde import gaussian_kde
from scipy.stats import norm


class KmerDistribution(object):
    def __init__(self, data_directory):
        self.data_directory = data_directory


class KmerHistogram(KmerDistribution):
    def __init__(self, data_directory, kmer):
        super(KmerHistogram, self).__init__(data_directory=data_directory)
        self.kmer = kmer
        self.data_file = data_directory + "{}_hist.txt".format(self.kmer)

        assert (os.path.isfile(self.data_file)), "No data for kmer {kmer}".format(kmer=self.kmer)

        self.histogram = None
        self.x_vals = None
        self.kde_pdf = None

        self.parse_histogram()

    def parse_histogram(self):
        self.histogram = np.loadtxt(self.data_file, dtype=np.float64)

    def parse_xvals(self, x_vals_file):
        self.x_vals = np.loadtxt(x_vals_file, dtype=np.float64)

    def make_kde(self, x_vals=None):
        kernel = gaussian_kde(self.histogram)
        ind = self.x_vals if self.x_vals is not None else x_vals
        assert ind is not None, "Need to provide x-values"
        self.kde_pdf = kernel.evaluate(ind)


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

