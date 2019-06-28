#!/usr/bin/env python

import sys
from Bio import SeqIO

for record in SeqIO.parse(sys.stdin, 'fasta'):
    print('>{id}\n{seq}'.format(id = record.id,
                                seq = record.seq[20:-20]))
