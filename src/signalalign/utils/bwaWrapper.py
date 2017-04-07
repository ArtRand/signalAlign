from __future__ import print_function

import sys
import os
import subprocess

import pysam
import numpy as np

from signalalign import _parseCigar


class GuideAlignment(object):
    def __init__(self, cigar=False, strand=False, reference_name=False):
        self.cigar = cigar
        self.strand = strand
        self.reference_name = reference_name

    def validate(self, reference_names):
        if self.cigar == "" or self.cigar is False or self.cigar is None:
            print("[GuideAlignment]Illegal cigar %s" % self.cigar)
            return False
        if self.strand not in ("+", "-"):
            print("[GuideAlignment]Illegal strand %s" % self.strand)
            return False
        if self.reference_name not in reference_names:
            print("[GuideAlignment]Reference %s not in reference names %s: " % (self.reference_name, reference_names))
            return False
        return True


class Bwa(object):
    """Wrapper of BWA aligner, requires bwa to be in path.
    Citation:
        Program: bwa (alignment via Burrows-Wheeler transformation)
        Contact: Heng Li <lh3@sanger.ac.uk>
    """
    def __init__(self, target):
        self.target = target
        self.db_handle = ''

    def build_index(self, destination, output=None):
        self.db_handle = destination + '/temp_bwaIndex'
        cmd = "bwa index -p {0} {1}".format(self.db_handle, self.target)
        if output is None:
            output = open(os.devnull, 'w')
        else:
            output = open(output, 'w')
        try:
            subprocess.check_call(cmd.split(), stdout=output, stderr=output)
            output.close()
            return True
        except subprocess.CalledProcessError:
            output.close()
            return False

    @staticmethod
    def suffixes():
        return [".amb", ".ann", ".bwt", ".pac", ".sa"]

    @staticmethod
    def align(bwa_index, query, output_sam_path, outerr=None):
        for suff in Bwa.suffixes():
            assert os.path.exists(bwa_index + suff),\
                "[Bwa:align] Didn't find index files {}".format(bwa_index + suff)
        assert os.path.exists(query), "[Bwa::align] Didn't find query file {}".format(query)
        cmd = "bwa mem -x ont2d {idx} {query}".format(idx=bwa_index, query=query)
        if outerr is None:
            outerr = open(os.devnull, 'w')
        else:
            outerr = open(outerr, 'w')
        try:
            with open(output_sam_path, 'w') as fH:
                fH.write(subprocess.check_output(cmd.split(), stderr=outerr))
            outerr.close()
            return True
        except subprocess.CalledProcessError:
            outerr.close()
            return False


def getBwaIndex(reference, dest, output=None):
    bwa = Bwa(reference)
    bwa.build_index(dest, output=output)
    bwa_ref_index = dest + "temp_bwaIndex"
    return bwa_ref_index


def generateGuideAlignment(bwa_index, query, temp_sam_path, target_regions=None):
    # type: (string, string, string, TargetRegions)
    """Aligns the read sequnece with BWA to get the guide alignment,
    returns the CIGAR (in exonerate format), the strand (plus or minus) and the
    contig mapped to if the read aligned. Returns (False, False, False) if there
    is a problem with any of the steps or if the read maps to a region not included
    within TargetRegions
    """
    class TargetRegions(object):
        def __init__(self, tsv, already_sorted=False):
            assert(os.stat(tsv).st_size != 0), "Empty regions file"

            self.region_array = np.loadtxt(tsv, usecols=(0, 1), dtype=np.int32)

            if len(self.region_array.shape) == 1:
                a = np.empty([1, 2], dtype=np.int32)
                a[0] = self.region_array
                self.region_array = a

            if not already_sorted:
                self.region_array = np.sort(self.region_array, axis=1)

        def check_aligned_region(self, left, right):
            if right < left:
                left, right = right, left
            for region in self.region_array:
                if (region[0] >= left) and (region[1] <= right):
                    return True
                else:
                    continue
            return False

    # align with bwa
    ok = Bwa.align(bwa_index=bwa_index, query=query, output_sam_path=temp_sam_path)
    if not ok:  # Bwa alignment fail
        return GuideAlignment()

    sam = pysam.Samfile(temp_sam_path, 'r')
    n_aligned_segments = 0
    query_name, flag, reference_name, reference_pos, sam_cigar = None, None, None, None, None

    for aligned_segment in sam:
        if not aligned_segment.is_secondary and not aligned_segment.is_unmapped:
            if n_aligned_segments == 0:
                query_name     = aligned_segment.qname
                flag           = aligned_segment.flag
                reference_name = sam.getrname(aligned_segment.rname)
                reference_pos  = aligned_segment.pos + 1  # pysam gives the 0-based leftmost start
                sam_cigar      = aligned_segment.cigarstring
            n_aligned_segments += 1

    if n_aligned_segments == 0:
        print("[exonerated_bwa_pysam]Read has no aligned segments")
        return GuideAlignment()

    if sam_cigar is None:
        print("[exonerated_bwa_pysam]DEBUG: query name: {qn} flag {fl} reference name {rn} "
              "reference pos {rp} sam cigar {cig} n_aligned {nal}"
              "".format(qn=query_name, fl=flag, rn=reference_name, rp=reference_pos, cig=sam_cigar,
                        nal=n_aligned_segments))

    if n_aligned_segments > 1:
        print("[exonerated_bwa_pysam]WARNING more than 1 mapping, taking the first one heuristically")

    try:
        query_start, query_end, reference_start, reference_end, cigar_string = _parseCigar(sam_cigar, reference_pos)
    except AssertionError, e:
        print("[generateGuideAlignment]ERROR %s" % e)
        return GuideAlignment()

    strand = ""
    assert(flag is not None), "[exonerated_bwa_pysam] ERROR flag is None"

    if int(flag) == 16:
        strand = "-"
        temp = reference_start
        reference_start = reference_end
        reference_end = temp
    if int(flag) == 0:
        strand = "+"
    elif int(flag) != 0 and int(flag) != 16:
        print("[exonerated_bwa_pysam]ERROR unexpected alignment flag {flag}, not continuing with signal alignment"
              " for {query}".format(flag=flag, query=query_name), file=sys.stderr)
        return GuideAlignment()

    assert(reference_name is not None), "[exonerated_bwa_pysam] ERROR reference_name is None"
    assert(query_name is not None), "[exonerated_bwa_pysam] ERROR query_name is None"

    completeCigarString = "cigar: %s %i %i + %s %i %i %s 1 %s" % (
        query_name, query_start, query_end, reference_name, reference_start, reference_end, strand, cigar_string)

    if target_regions is not None:
        target_regions = TargetRegions(target_regions)
        keep = target_regions.check_aligned_region(reference_start, reference_end)
        if keep is False:
            print("[exonerated_bwa_pysam]Read does not map witin the target regions, passing "
                  "on signal-level alignment", file=sys.stderr)
            return GuideAlignment()
        else:
            pass

    return GuideAlignment(completeCigarString, strand, reference_name)
