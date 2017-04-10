import os
import re


ALLOWED_FLAGS = (0, 16)

DEFAULT_TRAINMODELS_OPTIONS = {
    "fofn": None,
    "fast5_dir": None,
    "positions_file": None,
    "motif": None,
    "label": None,
}


def parseFofn(fofn_file):
    files = []
    with open(fofn_file, "r") as fH:
        for l in fH:
            files.append(l.strip())
    assert len(files) > 0, "parse_fofn: error, didn't find any files in file of files {}".format(fofn_file)
    return files


def _parseCigar(cigar_string, ref_start):
    assert(cigar_string is not None), "ERROR got cigar {}".format(cigar_string)
    assert(ref_start is not None)
    # use a regular expression to parse the string into operations and lengths
    cigar_tuples = re.findall(r'([0-9]+)([MIDNSHPX=])', cigar_string)

    clipping = {"S", "H"}
    alignment_operations = {"M", "I", "D"}

    # make some counters
    query_start = 0
    past_start = False
    query_end = 0
    reference_start = ref_start - 1  # fence posts adjustment
    reference_end = 0

    exonerated_cigar = " ".join(["%s %i" % (operation, int(length)) for length, operation in
                                 cigar_tuples if operation in alignment_operations])

    # this is how you calculate the reference map region
    for length, op in cigar_tuples:
        if op in clipping and past_start is False:
            query_start += int(length)
        if op == "M" or op == "D":
            reference_end += int(length)
            if past_start is False:
                past_start = True
        if op == "M" or op == "I":
            query_end += int(length)
            if past_start is False:
                past_start = True

    query_end = query_end + query_start
    reference_end = reference_end + reference_start

    return query_start, query_end, reference_start, reference_end, exonerated_cigar


def exonerateCigarWithStrandOrientation(aligned_segment, samfile):
    query_name = aligned_segment.qname
    flag       = aligned_segment.flag
    ref_name   = samfile.getrname(aligned_segment.rname)
    r_pos      = aligned_segment.pos + 1  # pysam does 0-based leftmost starts
    sam_cigar  = aligned_segment.cigarstring

    q_start, q_end, r_start, r_end, cigar_string = _parseCigar(sam_cigar, r_pos)

    if flag == 16:  # reverse mapping
        strand  = "-"
        tmp     = r_start
        r_start = r_end
        r_end   = tmp
    elif flag == 0:  # normal mapping
        strand = "+"
    else:  # not a primary alignment or some other u
        return None, False

    exonerate_cigar = "cigar: %s %i %i + %s %i %i %s 1 %s" % (query_name,
                                                              q_start,
                                                              q_end,
                                                              ref_name,
                                                              r_start,
                                                              r_end,
                                                              strand,
                                                              cigar_string)

    return exonerate_cigar, True


def defaultModelFromVersion(version, strand, pop1_complement=False):
    def default_template_model_from_version(version):
        supported_versions = ["1.15.0", "1.19.0", "1.20.0", "1.22.2", "1.22.4", "1.23.0"]
        assert version in supported_versions, "got version {}".format(version)
        version_index = supported_versions.index(version)
        if version_index <= 2:
            r7_3_default_template_model = "../models/testModelR73_acegot_template.model"
            assert os.path.exists(r7_3_default_template_model), "Didn't find default template R7.3 model"
            return r7_3_default_template_model
        elif version_index == 5:
            r94_default_template_model = "../models/testModelR9p4_acegt_template.model"
            assert os.path.exists(r94_default_template_model), "Didn't find default R9.4 model"
            return r94_default_template_model
        else:
            r9_default_template_model = "../models/testModelR9_template.model"
            assert os.path.exists(r9_default_template_model), "Didn't find default template R9 model"
            return r9_default_template_model

    def default_complement_model_from_version(version, pop1_complement=False):
        supported_versions = ["1.15.0", "1.19.0", "1.20.0", "1.22.2", "1.22.4"]
        assert version in supported_versions, "got version {}".format(version)
        version_index = supported_versions.index(version)

        if version_index <= 2:
            r7_3_default_complement_model = "../models/testModelR73_acegot_complement.model" if not pop1_complement \
                else "../models/testModelR9_complement_pop2.model"
            assert os.path.exists(r7_3_default_complement_model), "Didn't find default complement R7.3 model"
            return r7_3_default_complement_model
        else:
            r9_default_complement_model = "../models/testModelR9_complement.model"
            assert os.path.exists(r9_default_complement_model), "Didn't find default complement R9 model"
            return r9_default_complement_model

    if strand == "template":
        return default_template_model_from_version(version)
    elif strand == "complement":
        return default_complement_model_from_version(version, pop1_complement)
    else:
        return None
