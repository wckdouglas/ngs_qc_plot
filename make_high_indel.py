#!/usr/env python

import pandas as pd
import os
import sys
import numpy as np


if len(sys.argv) != 3:
    sys.exit('[usage] python %s <repeat_index table> <indel cutoff>')

ref_table = sys.argv[1]
indel_cut_off = int(sys.argv[2])


for gdf in pd.read_csv(ref_table, sep='\t', chunksize = 10000):
    for contig, contig_df in gdf.groupby('contig'): 
        df = contig_df\
            .assign(indel_index = lambda d: d.negative_index + d.positive_index) \
            .query('indel_index >= %i ' %indel_cut_off) 

        count = 0
        for i, base in df.iterrows():
            if base['negative_index'] == base['indel_index']:
                start = base['start']
                mononucleotide = base['fwd_base']
                indel_index = base['indel_index']
                taken_base = 1
            elif taken_base != indel_index and base['fwd_base'] == mononucleotide:
                taken_base += 1
            elif taken_base == indel_index:
                assert base['positive_index'] == indel_index and base['fwd_base'] == mononucleotide,'Wrong parsing'
                end = base['start']
                line = '{contig}\t{start}\t{end}\tIndel{id}\t{indel_index}\t+\t{mononucleotide}' \
                        .format(contig = base['contig'],
                                start = start, 
                                end = end, 
                                id = count,
                                indel_index = indel_index, 
                                mononucleotide = mononucleotide)
                print(line, file= sys.stdout)
                count += 1
            else:
                print(base)
