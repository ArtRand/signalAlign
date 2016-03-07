#!/usr/bin/env python
"""Write out the KL distance between two kmer models
"""
from __future__ import print_function
import os, sys
from vis_kmer_distributions import *
from scipy.stats import entropy
from itertools import product
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser (description=__doc__)

    parser.add_argument('--pk_dir', action='store', default=None, required=True, type=str, dest='pk_dir',
                        help="Path to experimental kmer distriutions")
    parser.add_argument('--out_file', action='store', default=None, required=True, type=str, dest='out',
                        help="place to put distances")
    args = parser.parse_args()
    return args


def main(args):
    args = parse_args()
    file_with_ont_model = "../../tests/minion_test_reads/C/" \
                          "makeson_PC_MA_286_R7.3_ZYMO_C_1_09_11_15_1714_1_ch1_file1_strand.fast5"
    assert os.path.exists(file_with_ont_model), "Didn't find ONT model containing file"
    assert os.path.exists(args.out) is not True, "Out file {} already exists".format(args.out)

    out_file = open(args.out, 'w')

    x_vals = np.linspace(30, 90, 600)

    kl_distances = []
    add_to_distances = kl_distances.append
    print("Collecting distances for {pk} against ONT table\n".format(pk=args.pk_dir), file=sys.stdout)
    for kmer in product("ACGT", repeat=6):
        kmer = ''.join(kmer)
        template_pdf, complement_pdf = plot_ont_distribution(kmer=kmer, fast5=file_with_ont_model, x_vals=x_vals)
        hdp_distribution = KmerHdpDistribution(data_directory=args.pk_dir, kmer=kmer)
        ent = entropy(pk=hdp_distribution.density, qk=template_pdf, base=2)
        add_to_distances(ent)
        print("finished with kmer {kmer} entropy {ent}".format(kmer=kmer, ent=ent), file=sys.stderr)
        out_file.write("{ent}\n".format(ent=ent))
    out_file.close()
    print("\nFinished collecting distances", file=sys.stdout)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
