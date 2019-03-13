#!/usr/bin/env python

import os
import sys
from collections import defaultdict
from goatools import obo_parser

if len(sys.argv) != 3:
    sys.exit('[usage] python %s <obo file> <gaf file>' %sys.argv[0])


go_terms = obo_parser.GODag(sys.argv[1])
pathways = defaultdict(set)
with open(sys.argv[2], 'r') as gaf:
    for line in gaf:
        if not line.startswith('!'):
            fields = line.strip().split('\t')
            GO_id = fields[4]
            gene_id = fields[1]
            pathways[GO_id].add(gene_id)

for go_id, genes in pathways.items():
    try:
        go_term = go_terms[go_id].name
    except KeyError:
        print(go_id)
        sys.exit()
    line = go_id + '\t' + go_term + '\t' + '\t'.join(list(genes))
    print(line)



