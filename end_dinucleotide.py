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
        .reset_index() \
        .fillna(0) \
        .query('base != "N"')

def plot_ends(df, figurename):
    df['first_base'] = map(lambda x: x[0] + '_', df['base'])
    df['second_base'] = map(lambda x: x[1] , df['base'])
    df = df[~df.base.str.contains('N')]
    with sns.plotting_context('paper',font_scale = 1.2), \
            sns.color_palette("husl", 8):
        p = sns.FacetGrid(data = df, col = 'read_end',
                    row = 'first_base',
                    hue = 'second_base', aspect = 1.5,
                    hue_order = list('ACTG'), 
                    row_order = ['A_','C_','T_','G_'])
    p.map(plt.plot, 'positions','base_fraction')
    p.add_legend(title = ' ')
    p.set_titles('{col_name}: {row_name}')
    p.set_axis_labels('Positions','Fraction')
    p.savefig(figurename)
    print 'Written %s ' %figurename
    return 0


def main():
    if len(sys.argv) != 3:
        sys.exit('[usage] python %s <bamfile> <outprefix>' %(sys.argv[0]))

    positions_consider = 20
    end_nucleotide_dict = defaultdict(lambda : defaultdict(lambda : defaultdict(int)))
    bam_file = sys.argv[1]
    outprefix = sys.argv[2]
    figurename = outprefix + '.pdf'
    tablename = outprefix + '.csv'
    positions = range(positions_consider)
    with pysam.Samfile(bam_file,'rb') as bam:
        for count, aln in enumerate(bam):
            condition_1 = (not aln.is_unmapped and not aln.is_supplementary)
            condition_2 = (not aln.is_duplicate and aln.mapping_quality > 1)
            if condition_1 and condition_2:
                sequence = str(aln.query_alignment_sequence)
                sequence = sequence if not aln.is_reverse else reverse_complement(sequence)
                read = "5'" if aln.is_read1 else "3'"
                sequence = sequence[:positions_consider] if aln.is_read1 else sequence[:positions_consider]
                sequence1 = sequence.translate(complement)[::-1] if aln.is_read2 else sequence
                sequence2 = sequence1[1:]
                for pos, base1, base2 in izip(positions, sequence1, sequence2):
                    end_nucleotide_dict[read][pos][base1 + base2] += 1
            if count % 10000000 == 0:
                print 'Parsed %i alignments' %(count)
    df = pd.concat([make_dataframe(end_nucleotide_dict, end) for end in ["5'","3'"]])
    plot_ends(df, figurename)
    df.to_csv(tablename, index=False)
    print 'Written %s' %tablename
    return 0



if __name__ == '__main__':
    main()
