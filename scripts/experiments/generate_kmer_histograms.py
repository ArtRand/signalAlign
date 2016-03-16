#!/usr/bin/env python

import os
import sys
from alignmentAnalysisLib import Kmer_histogram
from itertools import product
from argparse import ArgumentParser
from multiprocessing import Process, current_process, Manager


def parse_args():
    parser = ArgumentParser(description=__doc__)

    # query files
    parser.add_argument('--alignments', '-a', action='store',
                        dest='alns', required=False, type=str, default=None,
                        help="alignment files, add file extension")
    parser.add_argument('--number_of_assignments', '-n', action='store', type=int, default=10000,
                        dest='max_assignments',
                        help='total number of assignments to collect FOR EACH GROUP')
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=4, type=int, help="number of jobs to run concurrently")
    parser.add_argument('--threshold', '-t', action='store', type=float, default=0.25, dest='threshold')
    parser.add_argument('--out', '-o', action='store', type=str, required=True, dest='out')

    return parser.parse_args()


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
            k = Kmer_histogram(**f)
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
    """.format(alns=args.alns, maxNb=args.max_assignments, thresh=args.threshold)

    print start_message

    kmers_of_interest = []

    for kmer in product("ACTG", repeat=6):
        kmer = ''.join(kmer)
        if "C" in kmer:
            kmers_of_interest.append(kmer)
        else:
            continue

    workers = args.nb_jobs
    work_queue = Manager().Queue()
    done_queue = Manager().Queue()
    jobs = []

    template_directory = check_for_destination_directory(args.out, "template_hist/")
    complement_directory = check_for_destination_directory(args.out, "complement_hist/")
    if template_directory is None or complement_directory is None:
        print >> sys.stderr, "problem making destination directories"
        sys.exit(1)

    for strand, destination in zip(["t", "c"], [template_directory, complement_directory]):
        for kmer in kmers_of_interest:
            hist_args = {
                "path_to_alignments": args.alns,
                "kmer": kmer,
                "strand": strand,
                "threshold": args.threshold,
                "max_assignments": args.max_assignments,
                "out_dir": destination,
            }
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
