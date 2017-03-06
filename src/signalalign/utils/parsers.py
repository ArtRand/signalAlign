#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#Art Rand

"""
This is a module of FASTA/.qual and FASTQ parsers.

The functions in this module are intended for use with biological sequence data
(DNA, and amino acid) in the form of FASTA, .qual (which accompany FASTA),
and FASTQ files. There are two kinds of functions in this module, Readers and
Makers.

Readers:
In general, a reader function accepts an argument which points to the file(s)
to be parsed and then any options associated with that kind of file.

For example:
To parse a FASTQ file (ex, fake.fq64) that has phred-64/solexa quality encoding;

>>>read_fastq('fake.fq64', offset=64, solexa=True)

This function yields a tuple in with the internal format:
    ('sequence identifier', 'comment', 'sequence', [Q-value])

    'sequence identifier' = the string following the '@' or '>' in FASTA up
                            until the first space (' ') 
    'comment' = everything else on the first line after the sequence title

    'sequence' = the sequence (nucleotide, FASTA, FASTQ or amino acid FASTA)

    [Q-value] = all of the functions in this module convert to a common quality
                score which is defined as Q=-10(log10(Perror))

Makers:
These functions are intended to be used with the Reader functions included in 
this module. In general they accept an argument in the form of a generator 
(the Reader function), an output destination (default is sys.stdout) and any 
options associated with making that kind of file (ex. FASTQ offset).

For example:
To make a FASTA and .qual file from a FASTQ with phred-64 quality encoding that
prints to two files out.fasta and out.qual;

>>>make_fasta_and_qual(read_fastq('fake.fq64', offset=64), 
                       outfasta='out.fasta',
                       outqual='out.qual)

Limitations (and expansions):
None of these parsers check to make sure the files make sense. In other words,
if you have a FASTQ with letters that aren't in the alphabet associated with 
nucleotide sequence, the parser will not recognize them.  Similarly with quality
values, these parsers will not guess what encoding the input file is using or 
sanity check whether the quality values make sense.  

Basically, if you put garbage in, you'll get garbage out. But hopefully it will 
be correctly parsed and formatted garbage. 

Some edge conditions are accepted; no identifier string for example, and will 
simply propagate as empty strings for the title.  If you're only interested in 
converting from one format to another quickly, than this shouldn't be too much 
of a problem.  If there is no sequence data, these entries will propagate as 
empty entries. 
"""

import sys
import re
from math import log10
#from __future__ import print_function
from itertools import groupby
import argparse

def read_fasta(fasta_file, ignore_case=True):
    """
    Generic FASTA parser.  
    
    Input:
    Takes one argument, which can point to a file-type or a generator that
    yields one line at a time in FASTA format.  
    
    Output:
    An entry is identified by alternating '>sequence_id comment' then 'sequence'
    on a new line optionally followed by '>sequence2_id comment' and so forth.
    For each entry, a tuple ('title', 'comment', 'sequence') is yielded.  

    This function is not used directly in the 'maker' functions that follow, but
    is used in the read_fasta_with_quality() function which follows.     
    
    This code was adoped from 'brentp'
    URL: https://www.biostars.org/p/710/
    """ 
    
    file = open(fasta_file, 'r')
#    file=fasta_file
    
    # fasta_iter groups each line from the input fasta file by whether or not
    # the line starts with '>' , returning a tuple (bool, line) where bool
    # indicates if the line started with '>' and line is a string from the 
    # file.
    fasta_iter = (x[1] for x in groupby(file, 
                                        lambda line: line.startswith(">"))
                                        )
    # This for-loop takes the two groups of line(s) from fasta_iter. The ID
    # lines are seperated into the 'title' and the 'comment' The following group
    # of lines (the sequence) is joined into a single line and assigned to the
    # variable 'seq'. 'title', 'comment', and 'seq' are yielded as a tuple in
    # that order. 
    for lgroup in fasta_iter:
        lgroup = lgroup.next()[1:].strip()
        lgroup = lgroup.split(' ', 1)
        title = lgroup[0]
        comment = ''
        if len(lgroup) > 1:
            comment = lgroup[1]
        seq = "".join(s.strip() for s in fasta_iter.next())
        seq = seq.replace(" ", "")
        if ignore_case:
            seq = seq.upper()
        fasta = (title, comment, seq)
        yield fasta

def read_fasta_with_quality(fasta_file, qual_file):
    """
    Generic FASTA and .qual parser.  

    Input:
    Takes two arguments; fasta_file which can be a file-type or a generator that
    yields one line at a time in FASTA format and qual_file which contains the 
    quality data as whitespace separated integers for the associated FASTA file.

    Output:
    An entry is identified by alternating '>sequence_id comment' then 
    'sequence/qual-value' on a new line optionally followed by '>sequence2_id
    comment' and so forth.  For each entry a tuple ('title', 'comment',
    'sequence', [Q-value]) is yielded.  

    This function is used by the 'maker' functions that follow.    

    This code was adoped from brentp
    URL: https://www.biostars.org/p/710/
    """
    qual_file = open(qual_file, 'r')
    quals = {}
    # Identical to fast_iter in the read_fasta() function 
    qual_iter = (x[1] for x in groupby(qual_file, 
                                       lambda line: line.startswith(">"))
                                       )
    
    # This for-loop takes the two groups of lines from qual_iter. The ID lines
    # are seperated into the 'title' and the 'comment' and the 'comment' is 
    # discarded. The following group of lines (quality data) is put into a 
    # list of integers called 'q'. 'q' is assigned to the value in a dict with
    # 'title' as the key.
    for group in qual_iter:
        group = group.next()[1:].strip()
        group = group.split(' ', 1)
        qual_title = group[0]
        q = " ".join(s.strip() for s in qual_iter.next())
        q = [int(s) for s in q.split()]
        quals[qual_title] = q

    # The quality values (as a list) are appended to the tuple yeilded by the
    # read_fasta() function to give the internal format 
    # ('title', 'comment', 'seq', [Q-value]) and this tuple is yielded
    for fq in read_fasta(fasta_file):
        fq = fq + (quals.get(fq[0]), )
        yield fq

def read_fastq(fastq_file, offset=33, solexa=False): #add arguments
    """
    Generic FASTQ parser.  

    Input:
    Takes three arguments; fasta_file which can be a file-type or a generator
    that yields one line at a time in FASTQ format, an offset for phred-33 and 
    phred-64 (which is simply equal to whichever encoding the fastq_file has)
    and an optional solexa=True/False argument for phred-64/solexa encoding.

    Output:
    Each entry is identified by it's '@SEQID comment', 'sequence', 
    '+optional_id comment', 'quality characters'.  The quality characters are 
    converted into Q-values based on the offset/solexa arguments and the 
    'title', 'comment', 'sequence', and [Q-value] are yielded.

    This function is used by the 'maker' functions that follow.
    
    This code was adopted from 'Gareth Rees' on stackexchange
    URL": http://codereview.stackexchange.com/questions/32897/efficient-parsing-of-fastq
    """
    
    # These four lines are a set of regular expressions that are used to
    # identify the lines from the generator by whether or not they match
    # a pattern.
    at_seqname_re = re.compile(r'@(.+)$') 
    sequence_re = re.compile(r'[!-*,-~]*$')
    plus_seqname_re = re.compile(r'\+(.*)$')
    quality_re = re.compile(r'[!-ÿ]*$')

    lines = open(fastq_file)
    # This is the main loop of the function that iterates through each line of
    # the FASTQ file and groups the lines into 'title', 'comment', 'sequence'
    # and [Q-Value].
    for line in lines:
        seqname, comment = [], ''
        m = at_seqname_re.match(line)
        if not m: #make title empty string in the case of no title in file
            title = ''
        if m: #split indentifier line into 'title' and 'comment'
            seqname = m.group(1)
            seqname = seqname.split(' ', 1)
            title = seqname[0]
            comment = ''
        if len(seqname) > 1:
            comment = seqname[1]
        # identify lines that are sequence
        try:
            sequence = []
            for line in lines:
                m = sequence_re.match(line)
                if not m:
                    break
                sequence.append(m.group(0))
            if not sequence: 
                break
        # identify lines that are quality-identifier ('+SEQNAME')
            m = plus_seqname_re.match(line)
             # This is here to build the contingency that there is no '+'
#            if not m:
#                print "Expected +<seqname>"
#            if m.group(1) not in ['', title]:
#                print "Not the correct seqname for qual"
        # identify lines that are quality characters
            quality = []
            n = sum(map(len, sequence))
            while n > 0:
                line = next(lines)
                m = quality_re.match(line)
                if not m:
                    print >> sys.stderr,  "No quality characters"
                n -= len(m.group(0))
                if n < 0:
                    print >> sys.stderr, "quality longer than sequence"
                quality.append(m.group(0))
            qlist = list(''.join(quality))
            o = int(offset)
            qlist = [ord(i)-o for i in qlist] # this is where conversion from 
                                              # character to phred takes place
            if solexa == True: #this is where we convert from solexa to sanger
                qlist = [int(round(10*log10(10**(i/10.0) + 1))) for i in qlist]
            else:
                pass
            # this is left for debugging
#            print qlist
#            print seqname, ''.join(sequence), qlist
#            if not ''.join(sequence):
#                break
            fsq = (title, comment, ''.join(sequence), qlist)
#            print fsq
            yield fsq
            
        except StopIteration:
            print "End of input"

def make_fasta(input_fcn, outfasta=sys.stdout):
    """
    FASTA 'Maker' function.

    Input:
    Takes two arguments. The first, 'input_fcn' which a generator ('Read
    function') from above, although any generator that yields ('string',
    'string', 'string', [list of integers]) will work. The second, 'outfasta' 
    is the destination FASTA file or defaults to sys.stdout.

    Output:
    Default output is to sys.stdout.  Generates a string with the format:
    >'title' 'comment'
    sequence
    >'title2' 'comment2'
    ...
    
    NOTE: the quality data from the 'input_fcn' is discarded.  Also this 
    functions will not work directly with the read_fasta() function, so:
    >>>make_fasta(read_fasta(f.fa), output=newfasta.fa) will NOT work.
    A fasta_to_fasta function is in the works, the utility of which would only 
    be to re-format FASTA files for later use. 
    """
    # Check for non-default output
    if outfasta != sys.stdout:
        outfasta = open(outfasta, 'a')
    
    # For each tuple yielded by the input_fcn, write a FASTA-formatted string
    for title, comment, seq, qual in input_fcn:
        print >> outfasta, ">%s %s\n%s" % (title, comment, seq)

def make_fasta_and_qual(input_fcn, outfasta=sys.stdout, outqual=sys.stdout):
    """
    FASTA and .qual 'maker' function.

    Input:
    Takes three arguments. The first, 'input_fcn' which a generator ('Read
    function') from above, although any generator that yields ('string',
    'string', 'string', [list of integers]) will work. The second, 'outfasta'
    is the destination FASTA file or defaults to sys.stdout. The third,
    'outqual' is the destination quality file or defaults to sys.stdout 
    

    Output:
    Default output is to sys.stdout.  Generates a string with the format:
    >'title' 'comment'
    Sequence1
    >'title' 'comment'
    Q-values seperated by '  '
    >'title2', 'comment2'
    Sequence2
    ...
    
    NOTE: In the example both quality and sequence information are sent to 
    sys.stdout, but they can be routed to different output locations.
    """
    # Check for non-default output
    if outfasta != sys.stdout:
        outfasta = open(outfasta, 'a')
    if outqual != sys.stdout:
        outqual = open(outqual, 'a')

    # For each tuple yielded by the input_fcn, write a FASTA-formatted string
    # and a FASTA-formatted .qual file
    for title, comment, seq, qual in input_fcn:
        print >> outfasta, ">%s %s\n%s" % (title, comment, seq)
        print >> outqual, ">%s %s\n%s" % (title, comment,
                                          str(qual).strip('[]').
                                          replace(',', '  '))

def make_fastq(input_fcn, outfastq=sys.stdout, offset=33):
    """
    FASTQ 'maker' function.

    Input:
    Takes three arguments. The first, 'input_fcn' which a generator ('Read
    function') from above.  The second, 'outfastq' is the destination FASTQ file
    or defaults to sys.stdout. The third, 'offset' is the phred-33 or phred-64 
    offset.
    
    Output:
    Default output is to sys.stdout.  Generates a string with the format:
    @'title' 'comment'
    sequence
    +'title'
    Quality characters in phred-33 or phred-64 encoding
    @'title2' 'comment2'
    ...
    
    NOTE: there is no option to encode in phred-64/solexa. The 'title' string is
    added to the quality identifier line (+'title') even if it does not occur in
    the input file, the comment is not amended.  
    """
    # Check for non-default output
    if outfastq != sys.stdout:
        outfastq = open(outfastq, 'a')

    # For each tuple yielded by the input_fcn, write a FASTQ-formatted string
    for title, comment, seq, qual in input_fcn:
        qual = [chr(i+offset) for i in qual]
        print >> outfastq, "@%s %s\n%s\n+%s\n%s" % (title, comment, seq, 
                                                    title,
                                                  ''.join(qual))

