#!/usr/bin/env python2.7

import string
from collections import Counter

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
    trans_table = string.maketrans('ATGCatgc', 'TACGtacg')
    
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
        kmer = dna[i:(i+k)]
        if len(kmer) == k:
            kmer_count[kmer] += 1
    return kmer_count