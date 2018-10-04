#!/usr/bin/env python 

import pysam
import re
import sys


padded = 50

def change_header(header_line):
    seq_length = header_line.split(':')[-1]
    mod_length = str(int(seq_length) + padded * 2)
    return re.sub(seq_length + '$', mod_length, header_line)

def change_alignment(line):
    seq_start = line.split('\t')[3]
    mate_start = line.split('\t')[7]
    mod_start = str(int(seq_start) + padded)
    mate_mod_start = str(int(mate_start) + padded)
    line =  re.sub('\t' + seq_start + '\t', '\t' + mod_start + '\t', line)
    line =  re.sub('\t' + mate_start + '\t', '\t' + mate_mod_start + '\t', line)
    return line


def process_header(bam):
    for line in str(bam.header).split('\n'):
        if line.startswith('@SQ'):
            line = change_header(line)
        print(line)

with pysam.Samfile('-', 'r') as inbam:
    process_header(inbam)

    for aln in inbam:
        aln_string = aln.to_string()
        print(change_alignment(aln_string))

    



