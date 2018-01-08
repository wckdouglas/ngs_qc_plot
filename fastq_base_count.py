#!/usr/bin/env python

import sys
programname = sys.argv[0]
if len(sys.argv) < 2:
    sys.exit('python %s <file1> <file2> ...' %programname)
import matplotlib
matplotlib.use('Agg')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import Pool
import pyximport
pyximport.install()
from parse_fastq import read_file
sns.set_style('white')


def plot_seq(df):
	column_wrap_level = 2 if len(np.unique(df.sample_name)) > 1 else 1
	with sns.plotting_context('paper',font_scale=1.3):
		p = sns.FacetGrid(data = df, col = 'sample_name', hue='base',
						col_wrap = column_wrap_level, size=5,aspect = 2)
	p.map(plt.plot, 'position', 'norm_count')
	p.set_titles('{col_name}')
	p.set_xticklabels(rotation=30)
	p.add_legend()
	figurename = 'base_count_table.png'
	p.savefig(figurename)
	print 'Plotted: %s' %figurename


def norm_count(d):
    d['norm_count'] = np.true_divide(d['count'],np.sum(d['count'])) * 100
    return d

def main():
    files = sys.argv[1:]
    pool = Pool(24)
    tablename = 'base_count_table.tsv'
    dfs = pool.map(read_file, files)
    df = pd.concat(dfs, axis=0) \
		.to_csv(tablename,sep='\t',index=False)
    df = pd.read_table(tablename) \
		.groupby(['position','sample_name'], as_index=False) \
		.apply(norm_count)
    print 'Saved: %s' %tablename
    plot_seq(df)
    return 0

if __name__ == '__main__':
	main()
