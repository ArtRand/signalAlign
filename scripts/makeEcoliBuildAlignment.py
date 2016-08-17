#!/usr/bin/env python
""" Make supervised training alignment from E coli alignments
"""
import sys
import pandas as pd
from alignmentAnalysisLib import parse_alignment_file, cull_list_of_alignment_files
from signalAlignLib import parse_substitution_file
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser(description=__doc__)

    # query files
    parser.add_argument('--alignments', '-a', action='append', dest='alns', required=True, type=str,
                        help="alignment files to get 6-mer assignments from. Does not label assignments")
    parser.add_argument('--methylated', '-ma', action='append', dest='positive', required=True, type=str,
                        help="alignments to get labeled events from, will label 6-mers based on positions in "
                             "'methylated_positions file")
    parser.add_argument('--null', '-na', action='append', dest='null', required=False, type=str, default=None,
                        help="alignments to get null events from, same as '-ma' gets assignments based on "
                             "'null_positions' file, but does not label them. Optional.")
    parser.add_argument('--number_of_assignments', '-nb', action='store', type=int, default=10000,
                        dest='max_assignments', help='total number of canonical assignments to get, Default: 10000')
    parser.add_argument('--number_of_labeled_assignments', '-l', action='store', type=int, default=10000,
                        dest='max_labels', help='total number of assignments to collect for labeled and (optionally) '
                                                'null')
    parser.add_argument('--threshold', '-t', action='store', type=float, default=0.75, dest='threshold')
    parser.add_argument('--label', action='store', default="E", dest='label', help="Only collect assignments above this"
                                                                                   "posterior probability.")
    parser.add_argument('--methylated_positions', '-m', action='store', dest='positive_positions', required=True,
                        help="file with positions to label as 'label' ")
    parser.add_argument('--null_positions', '-n', action='store', dest='null_positions', required=False,
                        default=None, help="file with extra positions to add without labeling, and positions to"
                                           "ignore in the canonical alignments")
    parser.add_argument("--kmer_length", action='store', dest='kmer_length', required=True, type=int)
    parser.add_argument('--out', '-o', action='store', type=str, required=True, dest='out_file')

    return parser.parse_args()


def substitute_kmer(kmer, i, sub="E"):
    l = list(kmer)
    assert l[i] == "C"
    l[i] = sub
    return ''.join(l)


def get_reverse_complement_positions(positions):
    return [x - 3 for x in positions]


def get_hit_range(hit, kmer_length=6):
    return range(hit - (kmer_length - 1), hit + 1)[::-1]


def get_reverse_complement_hit_range(hit, kmer_length=6):
    return range(hit, hit + kmer_length)


def get_assignments(data, hit_range, strand, substitution, threshold):
    assignments = []
    for i, h in enumerate(hit_range):
        crit = data['ref_pos'].map(lambda x: x == h)
        selected = data[crit].ix[(data['strand'] == strand) & (data['posterior_prob'] >= threshold)]
        kmers = selected['kmer'].tolist()
        kmers = [substitute_kmer(kmer=x, i=i, sub=substitution) for x in kmers]
        assignment = pd.DataFrame({"kmer": kmers,
                                   "event_mean": selected['event_mean'],
                                   "posterior_prob": selected['posterior_prob'],
                                   "ref_pos": selected['ref_pos'],
                                   "event_index": selected['event_index'],
                                   "strand": strand})
        assignments.append(assignment)
    return pd.concat(assignments)


def get_assignment_table(data, hits, strand, get_range_function, substitution, threshold, kmer_length):
    assignment_table = []
    for hit in hits:
        hit_range = get_range_function(hit, kmer_length=kmer_length)
        assignments = get_assignments(data=data, hit_range=hit_range, strand=strand, substitution=substitution,
                                      threshold=threshold)
        assignment_table.append(assignments)
    return assignment_table


def get_labeled_assignmets_from_alignment(alignment, positions, forward, threshold, substitution, kmer_length):
    data = parse_alignment_file(alignment)
    rc_positions = [x - 3 for x in positions]

    crit = data['ref_pos'].map(lambda p: p in positions)
    rc_crit = data['ref_pos'].map(lambda p: p in rc_positions)

    hits = set(data[crit].ix[data['strand'] == 't']['ref_pos'])
    rc_hits = set(data[rc_crit].ix[data['strand'] == 't']['ref_pos'])

    if forward is True:
        assignments = get_assignment_table(data=data, hits=hits, strand='t', get_range_function=get_hit_range,
                                           substitution=substitution, threshold=threshold, kmer_length=kmer_length)
        rc_assignments = get_assignment_table(data=data, hits=rc_hits, strand='c',
                                              get_range_function=get_reverse_complement_hit_range,
                                              substitution=substitution, threshold=threshold, kmer_length=kmer_length)
    else:
        assignments = get_assignment_table(data=data, hits=rc_hits, strand='t',
                                           get_range_function=get_reverse_complement_hit_range,
                                           substitution=substitution, threshold=threshold, kmer_length=kmer_length)
        rc_assignments = get_assignment_table(data=data, hits=hits, strand='c', get_range_function=get_hit_range,
                                              substitution=substitution, threshold=threshold, kmer_length=kmer_length)

    return assignments + rc_assignments


def get_labeled_assignments(alignments, max_labels, positions, threshold, label, kmer_length=6):
    labeled_assignments = []  # bin for DataFrames of assignments
    total = 0  # count number of assignments we've gotten

    for alignment in alignments:
        if alignment.endswith(".backward.tsv"):
            forward = False
        else:
            forward = True
        # assignments is an array of DataFrames
        assignments = get_labeled_assignmets_from_alignment(alignment=alignment, positions=positions,
                                                            forward=forward, substitution=label,
                                                            threshold=threshold, kmer_length=kmer_length)
        assignments = [x for x in assignments if not x.empty]

        new_assignment_total = sum([x.shape[0] for x in assignments])
        print "Got {assignments} for {file}".format(assignments=new_assignment_total, file=alignment)

        total += sum([x.shape[0] for x in assignments])

        labeled_assignments += assignments
        if total >= max_labels:
            break

    if total < max_labels:
        print "Didn't get enough labels only got {notEnough} was aiming for {wanted} for label {label}" \
              "".format(notEnough=total, wanted=max_labels, label=label)
        return None

    return labeled_assignments


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
    # wrangle the positions we're calling positive and null
    pf, pb = parse_substitution_file(args.positive_positions)
    f_positive_substitution, forward_positive_sites = pf[0], pf[1]

    if args.null_positions is not None:
        nf, nb = parse_substitution_file(args.null_positions)
        f_null_substitution, forward_null_sites = nf[0], nf[1]
    else:
        forward_null_sites = []

    # randomly collect some alignments, and remove the events that correspond to positions we're going to label
    alignments = cull_list_of_alignment_files(args.alns)
    print "Getting canonical assignments"
    canonical_assignments = collect_assignments(alignments=alignments,
                                                threshold=args.threshold,
                                                max_assignments=args.max_assignments,
                                                positions=(forward_positive_sites + forward_null_sites))

    print "Got {} canonical assignments".format(canonical_assignments.shape[0])

    print "Getting labeled assignments"
    # get the list of alignments we're going to label
    positive_alignments = cull_list_of_alignment_files(args.positive)
    positive_assignments = get_labeled_assignments(alignments=positive_alignments, max_labels=args.max_labels,
                                                   positions=forward_positive_sites, threshold=args.threshold,
                                                   label=args.label, kmer_length=args.kmer_length)
    if positive_assignments is None:
        sys.exit(1)

    if args.null is not None:
        print "Getting null assignments"
        null_alignments = cull_list_of_alignment_files(args.null)
        null_assignments = get_labeled_assignments(alignments=null_alignments, max_labels=args.max_labels,
                                                   positions=forward_null_sites, threshold=args.threshold,
                                                   label="C")
        if null_assignments is None:
            sys.exit(1)

        total_labeled = sum([x.shape[0] for x in positive_assignments])
        total_null = sum([x.shape[0] for x in null_assignments])
        print "Got {labeled} labeled assignments and {null} null assignments".format(labeled=total_labeled,
                                                                                     null=total_null)
    else:
        null_assignments = None

    entry_line = "blank\t0\tblank\tblank\t{strand}\t0\t0.0\t0.0\t0.0\t{kmer}\t0.0\t0.0\t0.0\t{event}\t0.0\n"

    print "Writing to file {}".format(args.out_file)

    with open(args.out_file, 'w') as f:
        for row in canonical_assignments.itertuples():
            f.write(entry_line.format(strand=row[5], kmer=row[2], event=row[1]))

        for dataframe in positive_assignments:
            for row in dataframe.itertuples():
                f.write(entry_line.format(strand=row[6], kmer=row[3], event=row[2]))
        if null_assignments is not None:
            for dataframe in null_assignments:
                for row in dataframe.itertuples():
                    f.write(entry_line.format(strand=row[6], kmer=row[3], event=row[2]))

    print "DONE"

if __name__ == "__main__":
    sys.exit(main(sys.argv))