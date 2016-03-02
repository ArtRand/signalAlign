#!/usr/bin/env python
""" Make build alignments for starting a new HDP
"""
from __future__ import print_function
import glob
import os
import sys
import string
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from random import shuffle


def parse_args():
    parser = ArgumentParser(description=__doc__)

    # query files
    parser.add_argument('--C_alignments', '-C', action='store',
                        dest='C_alns', required=False, type=str, default=None,
                        help="C files")
    parser.add_argument('--mC_alignments', '-mC', action='store',
                        dest='mC_alns', required=False, type=str, default=None,
                        help="mC files")
    parser.add_argument('--hmC_alignments', '-hmC', action='store',
                        dest='hmC_alns', required=False, type=str, default=None,
                        help="hmC files")
    parser.add_argument('--number_of_assignments', '-n', action='store', type=int, default=10000,
                        dest='max_assignments',
                        help='total number of assignments to collect FOR EACH GROUP')
    parser.add_argument('--threshold', '-t', action='store', type=float, default=0.25, dest='threshold')
    parser.add_argument('--out', '-o', action='store', type=str, required=True, dest='out_file')

    return parser.parse_args()


def randomly_select_alignments(path_to_alignments):
    alignments = [x for x in glob.glob(path_to_alignments) if os.stat(x).st_size != 0]
    shuffle(alignments)
    return alignments


def collect_assignments(alignments, strand, threshold, max_assignments, transtable):
    if alignments is None:
        return None
    else:
        assignments_list = []
        add_to_assignments = assignments_list.append
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
                selected_rows = data.ix[(data['strand'] == strand) & (data['posterior_prob'] >= threshold)]
                total += selected_rows.shape[0]
                assignment_table = pd.DataFrame({"kmer": selected_rows['kmer'].str.translate(transtable),
                                                 "event_mean": selected_rows["event_mean"]})
                add_to_assignments(assignment_table)
            except:
                print("ERROR: problem with alignment {}".format(alignment))
                continue
            if total >= max_assignments:
                break
        assignments = pd.concat(assignments_list)
        return assignments


def make_build_alignment(c_alns, mc_alns, hmc_alns, strand, threshold, max_assignments):
    # translation tables for methylation
    C_trans_table = string.maketrans("C", "C")
    mC_trans_table = string.maketrans("C", "E")
    hmC_trans_table = string.maketrans("C", "O")

    C_table = collect_assignments(c_alns, strand, threshold, max_assignments, C_trans_table)
    mC_table = collect_assignments(mc_alns, strand, threshold, max_assignments, mC_trans_table)
    hmC_table = collect_assignments(hmc_alns, strand, threshold, max_assignments, hmC_trans_table)

    nb_c_assignments = C_table.shape[0] if C_table is not None else "None"
    nb_mc_assignments = mC_table.shape[0] if mC_table is not None else "None"
    nb_hmc_assignments = hmC_table.shape[0] if hmC_table is not None else "None"

    print("[buildAlignments] NOTICE: I found {C} C-assignments, {mC} mC-assignments, and {hmC} hmC-assignments "
          "for strand {strand}"
          "".format(C=nb_c_assignments, mC=nb_mc_assignments, hmC=nb_hmc_assignments, strand=strand),
          file=sys.stderr)
    tables = []

    for table in (C_table, mC_table, hmC_table):
        if table is None:
            continue
        else:
            tables.append(table)

    return pd.concat(tables)


def main(arguments):
    args = parse_args()

    C_alns = randomly_select_alignments(args.C_alns) if args.C_alns is not None else None
    mC_alns = randomly_select_alignments(args.mC_alns) if args.mC_alns is not None else None
    hmC_alns = randomly_select_alignments(args.hmC_alns) if args.hmC_alns is not None else None

    template_build_alignment = make_build_alignment(C_alns, mC_alns, hmC_alns, 't',
                                                    args.threshold, args.max_assignments)

    complement_build_alignment = make_build_alignment(C_alns, mC_alns, hmC_alns, 'c',
                                                      args.threshold, args.max_assignments)

    entry_line = "blank\t0\tblank\tblank\t{strand}\t0\t0.0\t0.0\t0.0\t{kmer}\t0.0\t0.0\t0.0\t{event}\t0.0\n"

    with open(args.out_file, 'w') as f:
        for row in template_build_alignment.itertuples():
            f.write(entry_line.format(strand="t", kmer=row[2], event=row[1]))
        for row in complement_build_alignment.itertuples():
            f.write(entry_line.format(strand="c", kmer=row[2], event=row[1]))

if __name__ == "__main__":
    sys.exit(main(sys.argv))