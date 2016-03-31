#!/usr/bin/env python
""" Make supervised training alignment from E coli alignments
"""
import sys
import glob
import os
import numpy as np
import pandas as pd
from random import shuffle
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser(description=__doc__)

    # query files
    parser.add_argument('--alignments', '-a', action='store', dest='alns', required=True, type=str,
                        help="alignment files")
    parser.add_argument('--labeled_alignments', '-la', action='append', dest='labeled', required=True, type=str,
                        help="alignments to get labeled events from")
    parser.add_argument('--number_of_assignments', '-n', action='store', type=int, default=10000,
                        dest='max_assignments',
                        help='total number of assignments to collect FOR EACH GROUP')
    parser.add_argument('--number_of_labeled_assignments', '-l', action='store', type=int, default=10000,
                        dest='max_labels', help='total number of assignments to collect FOR EACH GROUP')
    parser.add_argument('--threshold', '-t', action='store', type=float, default=0.75, dest='threshold')
    parser.add_argument('--positions', '-p', action='store', dest='positions', required=True)
    parser.add_argument('--out', '-o', action='store', type=str, required=True, dest='out_file')

    return parser.parse_args()


def parse_alignment_file(alignment_file):
    data = pd.read_table(alignment_file, usecols=(1, 4, 9, 12, 13),
                         dtype={'ref_pos': np.int64,
                                'strand': np.str,
                                'kmer': np.str,
                                'posterior_prob': np.float64,
                                'event_mean': np.float64},
                         header=None,
                         names=['ref_pos', 'strand', 'kmer', 'posterior_prob', 'event_mean'])
    return data


def parse_substitution_file(substitution_file):
    fH = open(substitution_file, 'r')
    line = fH.readline().split()
    forward_sub = line[0]
    forward_pos = map(np.int64, line[1:])
    line = fH.readline().split()
    backward_sub = line[0]
    backward_pos = map(np.int64, line[1:])
    return (forward_sub, forward_pos), (backward_sub, backward_pos)


def substitute_kmer(kmer, i, sub="E"):
    l = list(kmer)
    assert l[i] == "C"
    l[i] = sub
    return ''.join(l)

#def cull_assignments(alignment_df, positions, strand):
#    crit = data['ref_pos'].map(lambda x: x in positions)
#    cytosine_hits = set(data[crit].ix[data['strand'] == 't']['ref_pos'])


def get_reverse_complement_positions(positions):
    return [x - 3 for x in positions]


def get_hit_range(hit, kmer_length=6):
    return range(hit - (kmer_length - 1), hit + 1)[::-1]


def get_reverse_complement_hit_range(hit, kmer_length=6):
    return range(hit , hit + kmer_length)


def get_assignments(data, hit_range, strand, threshold=0.5):
    assignments = []
    for i, h in enumerate(hit_range):
        crit = data['ref_pos'].map(lambda x: x == h)
        selected = data[crit].ix[(data['strand'] == strand) & (data['posterior_prob'] >= threshold)]
        kmers = selected['kmer'].tolist()
        kmers = [substitute_kmer(x, i) for x in kmers]
        assignment = pd.DataFrame({"kmer": kmers,
                                   "event_mean": selected['event_mean'],
                                   "posterior_prob": selected['posterior_prob'],
                                   "ref_pos": selected['ref_pos'],
                                   "strand": strand})
        assignments.append(assignment)
    return pd.concat(assignments)


def get_assignment_table(data, hits, strand, get_range_function, threshold):
    assignment_table = []
    for hit in hits:
        hit_range = get_range_function(hit)
        assignments = get_assignments(data, hit_range, strand, threshold)
        assignment_table.append(assignments)
    return assignment_table


def get_labeled_assignmets_from_read(read, positions, forward, threshold, kmer_length=6):
    data = parse_alignment_file(read)
    rc_positions = [x - 3 for x in positions]

    crit = data['ref_pos'].map(lambda p: p in positions)
    rc_crit = data['ref_pos'].map(lambda p: p in rc_positions)

    hits = set(data[crit].ix[data['strand'] == 't']['ref_pos'])
    rc_hits = set(data[rc_crit].ix[data['strand'] == 't']['ref_pos'])

    if forward is True:
        assignments = get_assignment_table(data, hits, 't', get_hit_range, threshold)
        rc_assignments = get_assignment_table(data, rc_hits, 'c', get_reverse_complement_hit_range, threshold)
    else:
        assignments = get_assignment_table(data, rc_hits, 't', get_reverse_complement_hit_range, threshold)
        rc_assignments = get_assignment_table(data, hits, 'c', get_hit_range, threshold)

    return assignments + rc_assignments


def randomly_select_alignments(path_to_alignments):
    alignments = [x for x in glob.glob(path_to_alignments) if os.stat(x).st_size != 0]
    shuffle(alignments)
    return alignments


def collect_assignments(alignments, threshold, max_assignments, positions):
    if alignments is None:
        return None
    else:
        assignments_list = []
        add_to_assignments = assignments_list.append
        total = 0
        assert len(alignments) > 0, "Didn't find any alignments"
        for alignment in alignments:
            try:
                data = parse_alignment_file(alignment)
                crit = data['ref_pos'].map(lambda p: p not in positions)
                selected_rows = data[crit].ix[(data['posterior_prob'] >= threshold)]
                total += selected_rows.shape[0]
                assignment_table = pd.DataFrame({"kmer": selected_rows['kmer'],
                                                 "event_mean": selected_rows["event_mean"],
                                                 "ref_pos": selected_rows['ref_pos'],
                                                 "strand": selected_rows['strand'],
                                                 "posterior_prob": selected_rows['posterior_prob']})
                add_to_assignments(assignment_table)
            except:
                print("ERROR: problem with alignment {}".format(alignment))
                continue
            if total >= max_assignments:
                break
        assignments = pd.concat(assignments_list)
        return assignments


def main(args):
    args = parse_args()
    # randomly collect some alignments, and remove the events that correspont to methylated positions
    alignments = randomly_select_alignments(args.alns)
    f, b = parse_substitution_file(args.positions)
    forward_sub, forward_pos = f[0], f[1]
    print "Getting canonical assignments"
    canonical_assignments = collect_assignments(alignments,
                                                args.threshold,
                                                args.max_assignments,
                                                forward_pos)
    print "Got {} canonical assignments".format(canonical_assignments.shape[0])
    targeted_alignments = []
    for target_alignment in args.labeled:
        targeted_alignments += randomly_select_alignments(target_alignment)

    print "Getting labeled assignments"
    labeled_assignments = []
    total = 0
    while total < args.max_labels:
        for alignment in targeted_alignments:
            if alignment.endswith(".backward.tsv"):
                forward = False
            else:
                forward = True
            assignments = get_labeled_assignmets_from_read(alignment, forward_pos, forward, args.threshold)
            new_assignment_total = sum([x.shape[0] for x in assignments])
            print "Got {assignments} for {file}".format(assignments=new_assignment_total, file=alignment)
            total += sum([x.shape[0] for x in assignments])
            assignments = [x for x in assignments if not x.empty]
            labeled_assignments += assignments
        if total < args.max_labels:
            print "Didn't get enough labels only got {notEnough} was aiming for {wanted}, exiting" \
                  "".format(notEnough=total, wanted=args.max_labels)
            sys.exit(1)

    print "Got {} labeled assignments".format(total)

    entry_line = "blank\t0\tblank\tblank\t{strand}\t0\t0.0\t0.0\t0.0\t{kmer}\t0.0\t0.0\t0.0\t{event}\t0.0\n"

    print "Writing to file {}".format(args.out_file)

    with open(args.out_file, 'w') as f:
        for row in canonical_assignments.itertuples():
            f.write(entry_line.format(strand=row[5], kmer=row[2], event=row[1]))

        for frame in labeled_assignments:
            for row in frame.itertuples():
                f.write(entry_line.format(strand=row[5], kmer=row[2], event=row[1]))

    print "DONE"

if __name__ == "__main__":
    sys.exit(main(sys.argv))