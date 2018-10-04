#!/usr/bin/env python

from Bio import SeqIO
import sys

if len(sys.argv) != 2:
    sys.exit('[usage] python %s <fasta>' %sys.argv[0])

for record in SeqIO.parse(sys.argv[1],'fasta'):
    print('>{}\n{}'.format(record.id, 50*'N' + record.seq + 50*'N'))
