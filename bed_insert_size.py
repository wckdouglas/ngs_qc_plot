#!/bin/env python

import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np
import pysam
from pybedtools import BedTool
import seaborn as sns
import glob
import pandas as pd
import os
from multiprocessing import Pool
sns.set_style('white')

def parse_bed(bed_file):
    samplename = os.path.basename(bed_file).split('.')[0]
    print 'Extracting %s' %samplename
    fragSize = np.array([fragment.fields[4] for fragment in BedTool(bed_file)],dtype=np.int64)
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
    projectpath = '/scratch/02727/cdw2854/plasma_project'
    projectpath = '/scratch/02727/cdw2854/jurkatCells'
    datapath = projectpath + '/bedFiles'
    figurepath = datapath
    figurename = figurepath + '/insertSize.png'
    bed_files = glob.glob(datapath + '/*bed')
    dfs = Pool(24).map_async(parse_bed, bed_files).get()
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
