#!/usr/bin/env python

import glob
import os
import pandas as pd
import numpy as np
from random import shuffle


def randomly_select_alignments(path_to_alignments):
    alignments = [x for x in glob.glob(path_to_alignments) if os.stat(x).st_size != 0]
    shuffle(alignments)
    return alignments


class Kmer_histogram(object):
    def __init__(self, path_to_alignments, kmer, strand, threshold, max_assignments, out_dir):
        self.path_to_alignments = path_to_alignments
        self.kmer = kmer
        self.strand = strand
        self.threshold = threshold
        self.max_assignments = max_assignments
        self.output_directory = out_dir
        self.histogram = None
        self.n_points = 0

    def get_kmer_hist(self, alignments):
        kmer_hist = []
        add_to_hist = kmer_hist.append
        total = 0
        assert len(alignments) > 0, "Didn't find any alignments"
        for alignment in alignments:
            try:
                data = pd.read_table(alignment, usecols=(4, 9, 12, 13),
                                     dtype={'strand': np.str,
                                            'kmer': np.str,
                                            'posterior_prob': np.float64,
                                            'event_mean': np.float64},
                                     header=None,
                                     names=['strand', 'kmer', 'posterior_prob', 'event_mean'])

                selected_rows = data.ix[(data['strand'] == self.strand) &
                                        (data['posterior_prob'] >= self.threshold) &
                                        (data['kmer'] == self.kmer)]
                total += selected_rows.shape[0]
                aln_kmer_hist = pd.DataFrame({"event_mean": selected_rows["event_mean"]})
                add_to_hist(aln_kmer_hist)
            except Exception as e:
                print("ERROR: {exc} {aln}".format(exc=e, aln=alignment))
                continue
            if total >= self.max_assignments:
                    break
        self.histogram = pd.concat(kmer_hist)
        self.n_points = self.histogram.shape[0]
        return

    def output_kmer_hist(self):
        assert self.histogram is not None, "Histogram not initialized"
        out_file = self.output_directory + self.kmer + "_hist.txt"
        with open(out_file, 'w') as f:
            for row in self.histogram.itertuples():
                f.write("{mean}\n".format(mean=row[1]))
            f.write("\n")
            f.close()

    def run(self):
        alignments = randomly_select_alignments(self.path_to_alignments)
        self.get_kmer_hist(alignments)
        print "{n} points for {kmer} found".format(n=self.n_points, kmer=self.kmer)
        if self.n_points > 0:
            self.output_kmer_hist()
            return
        else:
            return

