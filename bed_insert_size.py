#!/bin/env python

import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np
import seaborn as sns
import glob
import pandas as pd
import os
from multiprocessing import Pool
sns.set_style('white')

def get_length(fragment):
    fields = fragment.split('\t')
    start, end = map(long, [fields[1], fields[2]])
    return abs(end - start)


def parse_bed(bed_file):
    samplename = os.path.basename(bed_file).split('.')[0]
    print 'Extracting %s' %samplename
    with open(bed_file,'r') as bed:
        fragSize = np.array([get_length(fragment) for fragment in bed],dtype=np.int64)
    fragSize = fragSize[fragSize<500]
    fragSize, count = np.unique(fragSize, return_counts=True)
    df = pd.DataFrame({'isize': fragSize, 'counts':count})
    df['samplename'] = samplename
    return df

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
    print 'plotted %s' %figurename
    return 0

def main():
    if len(sys.argv) != 3:
        sys.exit('[usage] python %s <threads> <bed_path>' %(sys.argv[0]))
    datapath = sys.argv[1]
    figurepath = datapath
    figurename = figurepath + '/insertSize.png'
    bed_files = glob.glob(datapath + '/*bed')
    dfs = Pool(int(sys.argv[2])).map_async(parse_bed, bed_files).get()
    df = pd.concat(dfs).\
            groupby(['samplename','isize'])\
            .sum()\
            .reset_index()\
            .groupby(['samplename'])\
            .apply(percetileDF)\
            .reset_index()
    plot(df, figurename)
    return 0

if __name__ == '__main__':
    main()
