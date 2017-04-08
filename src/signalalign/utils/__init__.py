import os

import string
from collections import Counter

import numpy as np
import pandas as pd

from signalalign.motif import getMotif
from signalalign.utils.parsers import read_fasta


def parseFofn(fofn_file):
    files = []
    with open(fofn_file, "r") as fH:
        for l in fH:
            files.append(l.strip())
    assert len(files) > 0, "parse_fofn: error, didn't find any files in file of files {}".format(fofn_file)
    return files


def reverse_complement(dna, reverse=True, complement=True):
    """
    Make the reverse complement of a DNA sequence. You can also just make the
    complement or the reverse strand (see options).

    Input: A DNA sequence containing 'ATGC' base pairs and wild card letters

    Output: DNA sequence as a string.

    Options: Specify reverse and/or complement as False to get the complement or
             reverse of the input DNA.  If both are False, input is returned.

    """

    # Make translation table
    trans_table = string.maketrans('ACGTMKRYBVDHNacgtmkrybvdhn',
                                   "TGCAKMYRVBHDNtgcakmyrvbhdn")
    # Make complement to DNA
    comp_dna = dna.translate(trans_table)

    # Output all as strings
    if reverse and complement:
        return comp_dna[::-1]
    if reverse and not complement:
        return dna[::-1]
    if complement and not reverse:
        return comp_dna
    if not complement and not reverse:
        return dna


def count_kmers(dna, k):
    """count the kmers of length k in a string"""
    kmer_count = Counter()
    for i in range(len(dna)):
        kmer = dna[i:(i + k)]
        if len(kmer) == k:
            kmer_count[kmer] += 1
    return kmer_count


class CustomAmbiguityPositions(object):
    def __init__(self, ambig_filepath):
        self.ambig_df = self.parseAmbiguityFile(ambig_filepath)

    @staticmethod
    def parseAmbiguityFile(ambig_filepath):
        """Parses a 'ambiguity position file' that should have the format:
            contig  position    strand  change_from change_to
        """
        return pd.read_table(ambig_filepath,
                             usecols=(0, 1, 2, 3, 4),
                             names=["contig", "position", "strand", "change_from", "change_to"],
                             dtype={"contig"      : np.str,
                                    "position"    : np.int,
                                    "strand"      : np.str,
                                    "change_from" : np.str,
                                    "change_to"   : np.str})

    def getForwardSequence(self, contig, raw_sequence):
        return self._get_substituted_sequence(contig, raw_sequence, "+")

    def getBackwardSequence(self, contig, raw_sequence):
        raw_sequence = reverse_complement(raw_sequence, reverse=False, complement=True)
        return self._get_substituted_sequence(contig, raw_sequence, "-")

    def _get_substituted_sequence(self, contig, raw_sequence, strand):
        contif_df = self._get_contig_positions(contig, strand)
        raw_sequence = list(raw_sequence)
        for _, row in contif_df.iterrows():
            if raw_sequence[row["position"]] != row["change_from"]:
                raise RuntimeError("[CustomAmbiguityPositions._get_substituted_sequence]Illegal substitution requesting "
                                   "change from %s to %s, row: %s" % (raw_sequence[row["position"]], row["change_to"], row))
            raw_sequence[row["position"]] = row["change_to"]
        return "".join(raw_sequence)

    def _get_contig_positions(self, contig, strand):
        return self.ambig_df.ix[(self.ambig_df["contig"] == contig) & (self.ambig_df["strand"] == strand)]


def processReferenceFasta(fasta, work_folder, motif_key=None, sub_char=None, positions_file=None):
    """loops over all of the contigs in the reference file, writes the forward and backward sequences
    as flat files (no headers or anything) for signalMachine, returns a dict that has the sequence
    names as keys and the paths to the processed sequence as keys
    """
    if positions_file is not None and motif_key is not None:
        raise RuntimeError("[processReferenceFasta]Cannot specify motif key and ambiguity position file")
    if positions_file is not None and sub_char is not None:
        raise RuntimeError("[processReferenceFasta]Cannot specify a substitution character and a ambiguity position file")

    if positions_file is not None:
        if not os.path.exists(positions_file):
            raise RuntimeError("[processReferenceFasta]Did not find ambiguity position file here: %s" %
                               positions_file)
        positions = CustomAmbiguityPositions(positions_file)
    else:
        positions = None

    ref_sequence_map = {}
    for header, comment, sequence in read_fasta(fasta):
        # the motif label allows us to make multiple copies of the reference with unique file names
        motif_lab = "" if motif_key is None else "%s." % motif_key
        # these are the paths to the flat files that have the references
        fw_path = work_folder.add_file_path("%s%s.%s.forward.txt" % (motif_lab, header, sub_char))
        bw_path = work_folder.add_file_path("%s%s.%s.backward.txt" % (motif_lab, header, sub_char))
        # signalAlign likes uppercase
        if motif_key is not None:
            motif, ok = getMotif(motif_key, sequence)
            if not ok:
                raise RuntimeError("[processReferenceFasta]Illegal motif key %s" % motif_key)
            fw_sequence = motif.forwardSubstitutedSequence(sub_char)
            bw_sequence = motif.complementSubstitutedSequence(sub_char)
        elif positions is not None:
            fw_sequence = positions.getForwardSequence(contig=header, raw_sequence=sequence.upper())
            bw_sequence = positions.getBackwardSequence(contig=header, raw_sequence=sequence.upper())
        else:
            fw_sequence = sequence.upper()
            bw_sequence = reverse_complement(fw_sequence, reverse=False, complement=True)

        with open(fw_path, 'w') as fH:
            fH.write("%s\n" % fw_sequence)
        with open(bw_path, 'w') as fH:
            fH.write("%s\n" % bw_sequence)

        ref_sequence_map[header] = {"forward": fw_path, "backward": bw_path}

    return ref_sequence_map
