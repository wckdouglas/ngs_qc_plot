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
    fields = line.split('\t')

    seq_start = fields[3]
    mate_start = fields[7]
    mod_start = str(int(seq_start) + padded)
    mate_mod_start = str(int(mate_start) + padded)
    line =  re.sub('\t' + seq_start + '\t', '\t' + mod_start + '\t', line)
    line =  re.sub('\t' + mate_start + '\t', '\t' + mate_mod_start + '\t', line)

    line = '\t'.join(fields[:3]) 
    line = line + '\t' + mod_start
    line = line + '\t' + '\t'.join(fields[4:7])
    line = line + '\t' + mate_mod_start 
    line = line + '\t' + '\t'.join(fields[8:])

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

    



