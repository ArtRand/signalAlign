from __future__ import print_function, division
import sys
import glob
import os
import pandas as pd
import numpy as np
from random import shuffle
from serviceCourse.parsers import read_fasta


def get_first_sequence(input_fasta):
    input_sequence = ""
    for header, comment, sequence in read_fasta(input_fasta):
        input_sequence += sequence
        break
    return input_sequence

def parse_alignment_file(alignment_file):
    data = pd.read_table(alignment_file, usecols=(1, 4, 5, 9, 12, 13),
                         dtype={'ref_pos': np.int64,
                                'strand': np.str,
                                'event_index': np.int64,
                                'kmer': np.str,
                                'posterior_prob': np.float64,
                                'event_mean': np.float64},
                         header=None,
                         names=['ref_pos', 'strand', 'event_index', 'kmer', 'posterior_prob', 'event_mean'])
    return data


def randomly_select_alignments(path_to_alignments):
    files = os.listdir(path_to_alignments)
    files = [f for f in files if f.endswith(".tsv")]
    files = [path_to_alignments + f for f in files]
    files = [f for f in files if os.path.isfile(f)]
    shuffle(files)
    print("[notice] sampling from {} files".format(len(files)), file=sys.stderr)
    return files


def get_alignments_from_directory(path_to_alignments):
    alignments = [x for x in glob.glob(path_to_alignments) if os.stat(x).st_size != 0]
    return alignments


def cull_list_of_alignment_files(list_of_directories):
    list_of_alignments = []
    for directory in list_of_directories:
        list_of_alignments += get_alignments_from_directory(directory)
    shuffle(list_of_alignments)
    return list_of_alignments


class KmerHistogram(object):
    def __init__(self, path_to_alignments, kmer, strand, threshold, max_assignments, ignore_positions, out_dir):
        self.path_to_alignments = path_to_alignments
        self.kmer = kmer
        self.strand = strand
        self.threshold = threshold
        self.max_assignments = max_assignments
        self.output_directory = out_dir
        self.ignore_positions = ignore_positions
        self.histogram = None
        self.n_points = 0

    def get_kmer_hist(self, alignments):
        kmer_hist = []
        add_to_hist = kmer_hist.append
        total = 0
        assert len(alignments) > 0, "Didn't find any alignments"
        for alignment in alignments:
            try:
                data = parse_alignment_file(alignment)
                if self.ignore_positions is not None:
                    crit = data['ref_pos'].map(lambda p: p not in self.ignore_positions)
                    selected_rows = data[crit].ix[(data['strand'] == self.strand) &
                                                  (data['posterior_prob'] >= self.threshold) &
                                                  (data['kmer'] == self.kmer)]
                else:
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
        alignments = cull_list_of_alignment_files(self.path_to_alignments)
        self.get_kmer_hist(alignments)
        print("{n} points for {kmer} found".format(n=self.n_points, kmer=self.kmer), file=sys.stdout)
        if self.n_points > 0:
            self.output_kmer_hist()
            return
        else:
            return


class CallMethylation(object):
    def __init__(self, sequence, alignment_file, degenerate_type, kmer_length, label=None, positions=None,
                 threshold=0.0, out_file=None):
        self.sequence = sequence
        self.forward = ".forward." in alignment_file
        self.alignment_file = alignment_file
        self.data = None
        self.probs = []
        self.template_calls = []
        self.complement_calls = []
        self.parse_alignment()
        self.out_file = out_file
        self.positions = positions
        self.degenerate = degenerate_type
        self.threshold = threshold
        self.kmer_length = kmer_length

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

    def get_range(self, position):
        return range(position - (self.kmer_length - 1), position + 1)

    def call_methyls(self, positions=None, threshold=0.0):
        if positions is None:
            template_sites = self.find_occurences("C") if self.forward is True else self.find_occurences("G")
            complement_sites = self.find_occurences("G") if self.forward is True else self.find_occurences("C")
        else:
            template_sites = positions['forward'] if self.forward is True else positions['backward']
            complement_sites = positions['backward'] if self.forward is True else positions['forward']

        def get_calls(sites, strand, regular_offset):
            for site in sites:
                # get the positions that an event can be aligned to and still report on this site
                positions = self.get_range(site)

                # select the rows in the dataFrame that have events aligned to this position
                crit = self.data['ref_index'].map(lambda x: x in positions)
                select = self.data[crit].ix[(self.data['strand'] == strand) & (self.data['prob'] >= threshold)]
                select = select.drop(select[['strand', 'ref_kmer']], axis=1)

                if select.empty:
                    continue

                marginal_probs = {"C": 0, "E": 0, "O": 0} if self.degenerate in [1, 2] \
                    else {"A": 0, "C": 0, "G": 0, "T": 0}

                for r in select.itertuples():
                    #kmer_length = len(r[5])
                    #offset = site - r[1] if regular_offset is True else 5 - (site - r[1])
                    offset = site - r[1] if regular_offset is True else (self.kmer_length - 1) - (site - r[1])
                    call = r[5][offset]
                    marginal_probs[call] += r[4]

                total_prob = sum(marginal_probs.values())

                for call in marginal_probs:
                    marginal_probs[call] /= total_prob

                self.probs.append((strand, site, marginal_probs))

        template_offset = True if self.forward is True else False
        complement_offset = False if self.forward is True else True
        get_calls(template_sites, 't', template_offset)
        get_calls(complement_sites, 'c', complement_offset)

    def write(self, out_file=None):
        self.call_methyls(positions=self.positions, threshold=self.threshold)
        fH = open(out_file, 'a') if self.out_file is None else open(self.out_file, 'a')

        file_name = self.alignment_file.split("/")[-1]

        def output_line():
            return "{site}\t{strand}\t{c}\t{mc}\t{hmc}\t{read}\n" if self.degenerate in [1, 2] \
                else "{site}\t{strand}\t{A}\t{C}\t{G}\t{T}\t{read}\n"

        line = output_line()

        for strand, site, prob in self.probs:
            if self.degenerate in [1, 2]:
                fH.write(line.format(site=site, strand=strand,
                                     c=prob["C"], mc=prob["E"], hmc=prob["O"],
                                     read=file_name))
            elif self.degenerate == 3:
                fH.write(line.format(site=site, strand=strand,
                                     A=prob["A"], C=prob["C"], G=prob["G"], T=prob["T"],
                                     read=file_name))
            else:
                sys.exit(1)
        return

    def run(self):
        self.write()
        return