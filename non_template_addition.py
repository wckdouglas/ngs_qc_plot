#!/usr/bin/env python

import pysam
import sys
import string
import numpy as np
import re
from collections import defaultdict
from operator import itemgetter

complementary = string.maketrans('ACTG','TGAC')
def reverse_complement(sequence):
    return sequence.translate(complementary)[::-1]

def cigar_to_seq(cigar):
    """
    Input a cigar string (eg. 1S16M5I10M4S)
    output a line: SMMMMMMMMMMMMMMMMIIIIIMMMMMMMMMMSSSS
    """
    cigarNum = np.array(re.findall('[0-9]+',cigar),dtype='int64')
    cigarStr = np.array(re.findall('[A-Z]',cigar),dtype='string')
    usable = np.in1d(cigarStr,np.array(['S','M','I'],dtype='string'))
    cigarStr = cigarStr[usable]
    cigarNum = cigarNum[usable]
    cigarSeq = ''
    for s, n in zip(cigarStr, cigarNum):
        cigarSeq += int(n)*str(s)
    return cigarSeq

def soft_clipped_seq(cigar,seq):
    soft_clipped = ''
    i = 0
    while i < len(cigar):
        c = cigar[i]
        if c == 'S':
            soft_clipped += seq[i]
        else:
            break
        i+=1
    return soft_clipped

def main():
    if len(sys.argv) != 2:
        sys.exit('[usage] python %s <bam_file>' %sys.argv[0])
    bam_file = sys.argv[1]
    clipped_pattern = defaultdict(int)
    mt_names = ['chrM','MT','Mt']
    with pysam.Samfile(bam_file,'rb') as in_bam:
        for aln in in_bam:
            if not aln.is_unmapped and 'S' in aln.cigarstring and aln.is_read1:
                chrom = in_bam.get_reference_name(aln.reference_id)
                if chrom not in mt_names:
                    cigar = cigar_to_seq(aln.cigarstring)
                    seq = aln.query_sequence
                    assert len(cigar) == len(seq), 'Wrong cigar'
                    if aln.is_reverse:
                        cigar = cigar[::-1]
                        seq = reverse_complement(seq)
                    if cigar.startswith('S'):
                        clipped = soft_clipped_seq(cigar, seq)
                        #print cigar + '\n' + clipped
                        clipped_pattern[clipped] += 1

    clipped_pattern = sorted(clipped_pattern.items(), 
                      key=itemgetter(1), reverse=True)
    for k, v in clipped_pattern:
        print '%s\t%i' %(k,v)


if __name__ == '__main__':
    main()
