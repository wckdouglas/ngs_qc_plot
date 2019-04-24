#!/usr/bin/env python

import pysam
from sequencing_tools.bam_tools import  split_cigar

def make_cigar(cs):
    cigar = split_cigar(cs)
    c = ''
    for n,b in zip(cigar[0], cigar[1]):
        if b != 'S':
            c = c + str(n) + str(b)
    return c



with pysam.Samfile('-') as inbam:
    with pysam.Samfile('-','wb',template = inbam) as outbam:
        for aln in inbam:
            if not aln.is_unmapped:
                cigar_str = aln.cigarstring
                seq = aln.query_alignment_sequence
                qual = aln.query_alignment_qualities
                aln.seq = seq
                aln.qual = ''.join(map(lambda x: chr(x + 33),  qual))
                aln.cigarstring = make_cigar(cigar_str)
                outbam.write(aln)
