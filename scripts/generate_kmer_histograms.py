#!/usr/bin/env python

import os
import sys
from alignmentAnalysisLib import KmerHistogram
from serviceCourse.parsers import read_fasta
from signalAlignLib import kmer_iterator, parse_substitution_file
from argparse import ArgumentParser
from multiprocessing import Process, current_process, Manager


def parse_args():
    parser = ArgumentParser(description=__doc__)

    # query files
    parser.add_argument('--alignments', '-a', action='append',
                        dest='alns', required=False, type=str, default=None,
                        help="alignment files, add file extension")
    parser.add_argument('--number_of_assignments', '-nb', action='store', type=int, default=100,
                        dest='max_assignments',
                        help='total number of points to collect')
    parser.add_argument('--ref', action='store', type=str, dest='ref', required=True,
                        help="Reference fasta file")
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=4, type=int, help="number of jobs to run concurrently")
    parser.add_argument('--threshold', '-t', action='store', type=float, default=0.8, dest='threshold')
    parser.add_argument('--ignore_positions', '-ig', action='store', dest='ignore', required=False, default=None,
                        help="file with positions to label as 'label' ")
    parser.add_argument('--out', '-o', action='store', type=str, required=True, dest='out')

    return parser.parse_args()


def get_sequence(path_to_fasta):
    seqs = []
    for header, comment, sequence in read_fasta(path_to_fasta):
        seqs.append(sequence)

    assert len(seqs) > 0, "ERROR parsing sequence {}".format(len(seqs))
    if len(seqs) > 1:
        print "Taking first sequence of {}".format(len(seqs))

    return seqs[0]


def check_for_destination_directory(working_directory_path, new_directory):
    try:
        if os.path.isdir(working_directory_path + new_directory):
            print "WARNING: destination for {} histograms already exists".format(new_directory)
        else:
            os.mkdir(working_directory_path + new_directory)
        return working_directory_path + new_directory
    except:
        return None


def histogram_runner(work_queue, done_queue):
    try:
        for f in iter(work_queue.get, 'STOP'):
            k = KmerHistogram(**f)
            k.run()
    except Exception, e:
        done_queue.put("%s failed with %s" % (current_process().name, e.message))


def main(args):
    args = parse_args()
    # get the kmers we want
    start_message = """
    # Generating Kmer histograms from alignments: {alns}
    # Getting a maximum of {maxNb} assignments
    # Using threshold: {thresh}
    # Ignoring events aligned to positions in {ignore}
    """.format(alns=args.alns, maxNb=args.max_assignments, thresh=args.threshold, ignore=args.ignore)

    print start_message

    kmers_of_interest = set()

    for kmer in kmer_iterator(get_sequence(args.ref), 6):
        kmers_of_interest.add(kmer)

    workers = args.nb_jobs
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()
    jobs = []

    template_directory = check_for_destination_directory(args.out, "template_hist/")
    complement_directory = check_for_destination_directory(args.out, "complement_hist/")
    if template_directory is None or complement_directory is None:
        print >> sys.stderr, "problem making destination directories"
        sys.exit(1)

    if args.ignore is not None:
        f, b = parse_substitution_file(args.ignore)
        ignore_positions = f[1]
    else:
        ignore_positions = None

    for strand, destination in zip(["t", "c"], [template_directory, complement_directory]):
        for kmer in kmers_of_interest:
            hist_args = {
                "path_to_alignments": args.alns,
                "kmer": kmer,
                "strand": strand,
                "threshold": args.threshold,
                "max_assignments": args.max_assignments,
                "out_dir": destination,
                "ignore_positions": ignore_positions,
            }
            #k = KmerHistogram(**hist_args)
            #k.run()
            work_queue.put(hist_args)

    for w in xrange(workers):
        p = Process(target=histogram_runner, args=(work_queue, done_queue))
        p.start()
        jobs.append(p)
        work_queue.put('STOP')

    for p in jobs:
        p.join()

    done_queue.put('STOP')

    print "\t# Finished Generating Kmer Histograms #"

if __name__ == "__main__":
    sys.exit(main(sys.argv))
