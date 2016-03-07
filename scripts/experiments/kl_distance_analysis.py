#!/usr/bin/env python
"""Write out the KL distance between two kmer models
"""
from __future__ import print_function
import os, sys
import numpy as np
from vis_kmer_distributions import *
from scipy.stats import entropy
from scipy.spatial.distance import euclidean
from itertools import product
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser (description=__doc__)

    parser.add_argument('--pk_dir', action='store', default=None, required=True, type=str, dest='pk_dir',
                        help="Path to experimental kmer distriutions")
    parser.add_argument('--out', action='store', default=None, required=True, type=str, dest='out',
                        help="place to put result files")
    args = parser.parse_args()
    return args


_SQRT2 = np.sqrt(2)


def hellinger2(p, q):
    return euclidean(np.sqrt(p), np.sqrt(q)) / _SQRT2


def main(args):
    args = parse_args()
    file_with_ont_model = "../../tests/minion_test_reads/C/" \
                          "makeson_PC_MA_286_R7.3_ZYMO_C_1_09_11_15_1714_1_ch1_file1_strand.fast5"
    assert os.path.exists(file_with_ont_model), "Didn't find ONT model containing file"
    kl_out_file_path = args.out + "kl_distance.txt"
    hd_out_file_path = args.out + "hellinger_distance.txt"
    assert os.path.exists(kl_out_file_path) is not True, "Out file {} already exists".format(kl_out_file_path)
    assert os.path.exists(hd_out_file_path) is not True, "Out file {} already exists".format(hd_out_file_path)

    kl_out = open(kl_out_file_path, 'w')
    hd_out = open(hd_out_file_path, 'w')

    x_vals = np.linspace(30, 90, 600)

    print("Collecting distances for {pk} against ONT table\n".format(pk=args.pk_dir), file=sys.stdout)
    for kmer in product("ACGT", repeat=6):
        kmer = ''.join(kmer)
        template_pdf, complement_pdf = plot_ont_distribution(kmer=kmer, fast5=file_with_ont_model, x_vals=x_vals)
        hdp_distribution = KmerHdpDistribution(data_directory=args.pk_dir, kmer=kmer)
        ent = entropy(pk=hdp_distribution.density, qk=template_pdf, base=2)
        h_distance = hellinger2(p=hdp_distribution.density, q=template_pdf)
        print("finished with kmer {kmer} entropy {ent} hellinger distance {hd}"
              "".format(kmer=kmer, ent=ent, hd=h_distance), file=sys.stderr)
        kl_out.write("{ent}\n".format(ent=ent))
        hd_out.write("{hd}\n".format(hd=h_distance))

    kl_out.close()
    hd_out.close()
    print("\nFinished collecting distances", file=sys.stdout)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
