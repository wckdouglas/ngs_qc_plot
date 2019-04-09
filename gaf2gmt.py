#!/usr/bin/env python

import os
import sys
from collections import defaultdict
from goatools import obo_parser

def make_go_term_dict(obo):
    go_terms = obo_parser.OBOReader(sys.argv[1])
    go_term_dict = {}
    for go in go_terms:
        go_term_dict[go.id] = go.name
        if go.alt_ids:
            for id in go.alt_ids:
                go_term_dict[id] = go.name
    return go_term_dict


if len(sys.argv) != 4:
    sys.exit('[usage] python %s <obo file> <gaf file> <gmt file>' %sys.argv[0])


go_term_dict = make_go_term_dict(sys.argv[1])
pathways = defaultdict(set)
with open(sys.argv[2], 'r') as gaf:
    for line in gaf:
        if not line.startswith('!'):
            fields = line.strip().split('\t')
            GO_id = fields[4]
            gene_id = fields[1]
            pathways[GO_id].add(gene_id)

with open(sys.argv[3], 'w') as gmt:
    for go_id, genes in pathways.items():
        go_term = go_term_dict[go_id]
        line = go_term + '\t' + go_id + '\t' + '\t'.join(list(genes))
        print(line, file = gmt)



