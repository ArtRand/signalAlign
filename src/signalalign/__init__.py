import re

ALLOWED_FLAGS = (0, 16)


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
