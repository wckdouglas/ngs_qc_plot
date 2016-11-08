#!/bin/env python

from matplotlib import use
use('Agg')
import pylab as plt
import numpy as np
import pysam
import seaborn as sns
import glob
import pandas as pd
import os
from multiprocessing import Pool
sns.set_style('white')

def filter_align(aln):
    if aln.is_proper_pair and not aln.is_secondary \
            and aln.is_read1 and not aln.is_supplementary:
        return abs(aln.isize)
    else:
        return 0

def parseBam(bamFile):
    samplename = os.path.basename(bamFile).split('.')[0]
    print 'Extracting %s' %samplename
    with pysam.Samfile(bamFile, 'rb') as bam:
        fragSize = np.abs([filter_align(aln) for aln in bam])
    fragSize = fragSize[fragSize>20]
    fragSize, count = np.unique(fragSize, return_counts=True)
    df = pd.DataFrame({'isize': fragSize, 'counts':count})
    df['samplename'] = samplename
    return df

def percetileDF(d):
    d['counts'] = np.true_divide(d['counts'],np.max(d['counts']))
    return d

def plot(df, figurename):
    with sns.plotting_context('paper',font_scale=1.4):
        p = sns.FacetGrid(data = df, col = 'samplename', col_wrap = 4, sharey=False)
    p.map(plt.bar,'isize','counts')
    p.set_titles('{col_name}')
    p.set_xticklabels(rotation=60)
	p.fig.text(x=0, y=0.7, s='Normalized Count', rotation=90)
    p.fig.text(x=0.5, y = 0, s='Fragment Size (nt)')
    p.savefig(figurename)
    print 'plotted %s' %figurename
    return 0

def main():
    if len(sys.argv) != 2:
        sys.exit('[usage] python %s <bam_path>' %(sys.argv[0]))
    datapath = sys.argv[1]
    figurepath = datapath + '/figures'
    if not os.path.isdir(figurepath):
        os.mkdir(figurepath)
    figurename = figurepath + '/insertSize.png'
    tablename = figurename.replace('.png','.tsv')
    bamFiles = glob.glob(datapath + '/P*.bam')
    p = Pool(24)
    dfs = p.map(parseBam, bamFiles)
    df = pd.concat(dfs)\
        .groupby(['samplename','isize'])\
	.agg({'counts':np.sum})\
        .reset_index()\
        .groupby(['samplename'])\
        .apply(percetileDF)\
        .reset_index()
    p.close()
    p.join()
    df.to_csv(tablename,sep='\t',index=False)
    plot(df, figurename)
    return 0


if __name__ == '__main__':
    main()
