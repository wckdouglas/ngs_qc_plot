#!/usr/bin/env python

import sys
programname = sys.argv[0]
if len(sys.argv) < 2:
    sys.exit('python %s <file1> <file2> ...' %programname)
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import matplotlib
matplotlib.use('Agg')
import gzip
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import Pool
sns.set_style('white')

def seq_base_dict(n_bases):
	seq_dict = {base: np.zeros(n_bases) for base in list('ACTGN')}
	return seq_dict

def parse_seq(sequence, seq_dict):
	for i, base in enumerate(sequence):
		seq_dict[base][i] += 1
	return seq_dict

def open_file(filename):
	if filename.endswith('gz') or filename.endswith('gzip'):
		return gzip.open(filename,'r')
	elif filename.endswith('fastq') or filename.endswith('fq'):
		return open(filename,'r')
	else:
		sys.exit('Unrecognized file type!')

def read_file(args):
    filename, n_bases = args
    print 'Running %s' %filename
    samplename = os.path.basename(filename).split('.')[0]
    seq_dict = seq_base_dict(n_bases)
    with open_file(filename) as fastq:
    	i = 0
    	for name, seq, qual in FastqGeneralIterator(fastq):
        	seq_dict = parse_seq(seq, seq_dict)
		i += 1
		if i % 100000 == 0:
	    	    print 'Parsed: %i records in %s' %(i, filename)
    df = pd.DataFrame(seq_dict) \
		.assign(sample_name = samplename) \
		.assign(position = np.arange(n_bases) + 1)
    return df

def plot_seq(df):
	column_wrap_level = 2 if len(np.unique(df.sample_name)) > 1 else 1
	with sns.plotting_context('paper',font_scale=1.3):
		p = sns.FacetGrid(data = df, col = 'sample_name', hue='base',
						col_wrap = column_wrap_level, size=5,aspect = 2)
	p.map(plt.plot, 'position', 'count')
	p.set_titles('{col_name}')
	p.set_xticklabels(rotation=30)
	p.add_legend()
	figurename = 'base_count_table.png'
	p.savefig(figurename)
	print 'Plotted: %s' %figurename

def norm_count(df):
	df['count'] = np.true_divide(df['count'],np.sum(df['count'])) * 100
	return df

def main():
	files = sys.argv[1:]
	n_bases = 200
	pool = Pool(24)
	tablename = 'base_count_table.tsv'
	dfs = pool.map(read_file, [(filename, n_bases) for filename in files])
	pd.concat(dfs, axis=0)\
		.pipe(pd.melt, id_vars=['position','sample_name'], value_vars=list('ACTGN'),
						var_name='base',value_name = 'count') \
		.groupby(['position','sample_name'])\
		.apply(norm_count)\
		.reset_index()\
		.to_csv(tablename,sep='\t',index=False)
	df = pd.read_table(tablename, sep='\t')
	print 'Saved: %s' %tablename
	plot_seq(df)
	return 0

if __name__ == '__main__':
	main()
