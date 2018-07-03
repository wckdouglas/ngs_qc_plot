#!/bin/env python

from __future__ import print_function
from matplotlib import use as mpl_use
mpl_use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import glob
import pandas as pd
import os
import sys
from multiprocessing import Pool
import pyximport
pyximport.install()
from parse_bed import parse_bed
sns.set_style('white')

def get_length(fragment):
    fields = fragment.split('\t')
    start, end = map(long, [fields[1], fields[2]])
    return abs(end - start)


def percetileDF(df):
    df['counts'] = np.true_divide(df['counts'],np.sum(df['counts'])) * 100
    return df

def plot(df, figurename):
    with sns.plotting_context('paper',font_scale=1.4):
        p = sns.FacetGrid(data = df, col = 'samplename', col_wrap = 2, sharey=False)
        p.map(plt.bar,'isize','counts')
        p.set(xlabel='Fragment Size (nt)')
        p.set(ylabel='Percentage of Reads')
        p.set_titles('{col_name}')
        p.savefig(figurename)
    print('plotted %s' %figurename)
    return 0

def main():
    if len(sys.argv) != 3:
        sys.exit('[usage] python %s <bed_path> <threads>' %(sys.argv[0]))
    datapath = sys.argv[1]
    figurepath = datapath
    figurename = figurepath + '/insertSize.png'
    bed_files = glob.glob(datapath + '/*bed')
    bed_files.extend(glob.glob(datapath + '/*bed.gz'))
    p = Pool(int(sys.argv[2]))
    dfs = p.map(parse_bed, bed_files)
    p.close()
    p.join()
    df = pd.concat(map(pd.read_csv, dfs)).\
            groupby(['samplename','isize'])\
            .sum()\
            .reset_index()\
            .groupby(['samplename'])\
            .apply(percetileDF)\
            .reset_index()
    tablename = figurename.replace('.png','.tsv')
    df.to_csv(tablename, index=False, sep='\t')
    print('Saved %s' %(tablename))
    plot(df, figurename)
    [os.remove(f) for f in dfs]
    return 0

if __name__ == '__main__':
    main()
