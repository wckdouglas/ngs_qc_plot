#!/usr/bin/env python

from Bio import SeqIO
import sys

if len(sys.argv) > 2:
    sys.exit('[usage] python %s <fasta>' %sys.argv[0])

elif len(sys.argv) == 1:
    infile = sys.stdin

else:
    infile = sys.argv[1]

inseq = SeqIO.parse(infile, 'fasta')


for record in inseq:
    print('>{}\n{}'.format(record.id, 50*'N' + record.seq + 50*'N'))
