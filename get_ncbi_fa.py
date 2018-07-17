#!/usr/bin/env python

from Bio import Entrez, SeqIO
import sys
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

if len(sys.argv) < 2:
    sys.exit('[usage] python %s <ncbi_id1> <ncbi_id2> ...')

IDs= sys.argv[1:]

Entrez.email = 'wckdouglas@gmail.com'

for id in IDs:
    id = id.split(' ')[0]
    handle = Entrez.efetch(db="nucleotide", id=id, 
                           rettype="fasta", retmode="text")
    record = handle.read()
    print ('>{id}\n{seq}'.format(id = id, seq='\n'.join(record.split('\n')[1:])))

