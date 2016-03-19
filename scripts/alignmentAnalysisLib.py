#!/usr/bin/env python

from __future__ import print_function, division
import sys
import glob
import os
import pandas as pd
import numpy as np
from random import shuffle


def randomly_select_alignments(path_to_alignments):
    files = os.listdir(path_to_alignments)
    files = [f for f in files if f.endswith(".tsv")]
    files = [path_to_alignments + f for f in files]
    files = [f for f in files if os.path.isfile(f)]
    shuffle(files)
    print("[notice] sampling from {} files".format(len(files)), file=sys.stderr)
    return files


class KmerHistogram(object):
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
        print("{n} points for {kmer} found".format(n=self.n_points, kmer=self.kmer), file=sys.stdout)
        if self.n_points > 0:
            self.output_kmer_hist()
            return
        else:
            return


class CallMethylation(object):
    def __init__(self, sequence, alignment_file, forward, label=None, test=False):
        self.sequence = sequence
        self.forward = forward
        self.alignment_file = alignment_file
        self.data = None
        self.template_calls = []
        self.complement_calls = []
        self.parse_alignment()
        if label is not None:
            self.label = label

    def parse_alignment(self):
        # todo make this try, except
        self.data = pd.read_table(self.alignment_file,
                                  usecols=(1, 2, 4, 5, 6, 7, 8),
                                  header=None,
                                  names=['ref_index', 'ref_kmer', 'strand', 'event_index',
                                         'match_kmer', 'prob', 'path_kmer'],
                                  dtype={'ref_index': np.int64,
                                         'ref_kmer': np.str,
                                         'strand': np.str,
                                         'event_index': np.int64,
                                         'match_kmer': np.str,
                                         'prob': np.float64,
                                         'path_kmer': np.str})
        assert self.data is not None, "problem parsing alignment file {}".format(self.alignment_file)

    def find_occurences(self, ch):
        return [i for i, letter in enumerate(self.sequence) if letter == ch]

    @staticmethod
    def get_range(position):
        return range(position - 5, position + 1)

    def call_methyls(self):
        template_sites = self.find_occurences("C") if self.forward is True else self.find_occurences("G")
        complement_sites = self.find_occurences("G") if self.forward is True else self.find_occurences("C")

        def get_calls(sites, call_bin, strand, regular_offset):
            for site in sites:
                # get the positions that an event can be aligned to and still report of this site
                positions = self.get_range(site)

                # select the rows in the dataFrame that have events aligned to this position
                crit = self.data['ref_index'].map(lambda x: x in positions)
                select = self.data[crit].ix[self.data['strand'] == strand]
                select = select.drop(select[['strand', 'ref_kmer']], axis=1)

                if select.empty:
                    continue

                probs = {
                    "C": 0,
                    "E": 0,
                    "O": 0,
                }

                for r in select.itertuples():
                    offset = site - r[1] if regular_offset is True else 5 - (site - r[1])
                    call = r[5][offset]
                    probs[call] += r[4]

                total_prob = sum(probs.values())

                for call in probs:
                    probs[call] /= total_prob

                call_bin.append(max(probs, key=probs.get))

        template_offset = True if self.forward is True else False
        complement_offset = False if self.forward is True else True
        get_calls(template_sites, self.template_calls, 't', template_offset)
        get_calls(complement_sites, self.complement_calls, 'c', complement_offset)

    def test(self):
        self.call_methyls()
        correct_calls = {
            "template": self.template_calls.count(self.label),
            "complement": self.complement_calls.count(self.label)
        }
        template_acc = self.template_calls.count(self.label) / \
                       len(self.template_calls) if len(self.template_calls) > 0 else 0
        complement_acc = self.complement_calls.count(self.label) / \
                         len(self.complement_calls) if len(self.complement_calls) > 0 else 0
        file_name = self.alignment_file.split("/")[-1]
        print("{file}\t{t}\t{c}".format(t=template_acc, c=complement_acc, file=file_name), file=sys.stderr)
        return correct_calls
