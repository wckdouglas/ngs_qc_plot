#!/usr/bin/env python

from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
import pysam
import numpy as np
from collections import defaultdict
import seaborn as sns
import pandas as pd
import string
from itertools import izip
import sys
import re

complement = string.maketrans('ACTGN','TGACN')
def reverse_complement(seq):
    return seq.translate(complement)[::-1]

def norm_data(d):
    d['base_fraction'] = np.true_divide(d['base_count'], d.base_count.sum())
    return d

def make_dataframe(nucleotide_dict, end):
    return pd.DataFrame.from_dict(nucleotide_dict[end], orient = 'index') \
        .reset_index() \
        .rename(columns = {'index':'positions'})\
        .assign(read_end = end)  \
        .pipe(pd.melt, id_vars = ['positions','read_end'],
                value_name = 'base_count', var_name='base')\
        .groupby(['read_end','positions']) \
        .apply(norm_data) \
        .reset_index()\
        .fillna(0) \
        .query('base != "N"')\
        .drop('index',axis=1)

def plot_ends(df, figurename):
    with sns.plotting_context('paper',font_scale = 1.2), \
            sns.color_palette("husl", 8):
        p = sns.FacetGrid(data = df, col = 'read_end',
                      hue = 'base', aspect = 1.5)
    p.map(plt.plot, 'positions','base_fraction')
    p.add_legend()
    p.set_titles('{col_name}')
    p.set_axis_labels('Positions','Fraction')
    p.savefig(figurename, transparent=True)
    print 'Written %s ' %figurename
    return 0


def good_cigar(cigar):
    cigar = str(cigar)
    return re.findall('[MHSID]',cigar) == ['M']


def extract_nucleotides(bam, positions_consider):
    end_nucleotide_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    positions = range(positions_consider)
    for count, aln in enumerate(bam):
        condition_1 = (not aln.is_unmapped and not aln.is_supplementary)
        condition_2 = (not aln.is_duplicate and aln.mapping_quality > 1)
        condition_3 = good_cigar(aln.cigarstring)
        if condition_1 and condition_2 and condition_3:
            #sequence = str(aln.query_alignment_sequence)
            read = "5'" if aln.is_read1 else "3'"
            sequence = str(aln.query_sequence) if read == "5'" else str(aln.query_sequence)
            sequence = sequence[:positions_consider] if aln.is_read1 else reverse_complement(sequence[:positions_consider])
            for pos, base in izip(positions, sequence):
                end_nucleotide_dict[read][pos][base] += 1
        if count % 10000000 == 0:
            print 'Parsed %i alignments' %(count)
    return end_nucleotide_dict


def main():
    if len(sys.argv) != 3:
        sys.exit('[usage] python %s <bamfile> <outprefix>' %(sys.argv[0]))

    positions_consider = 20
    bam_file = sys.argv[1]
    outprefix = sys.argv[2]
    figurename = outprefix + '.pdf'
    tablename = outprefix + '.csv'
    with pysam.Samfile(bam_file,'rb') as bam:
        end_nucleotide_dict = extract_nucleotides(bam, positions_consider)
    df = pd.concat([make_dataframe(end_nucleotide_dict, end) for end in ["5'","3'"]])
    plot_ends(df, figurename)
    df.to_csv(tablename, index=False)
    print 'Written %s' %tablename
    return 0


if __name__ == '__main__':
    main()
