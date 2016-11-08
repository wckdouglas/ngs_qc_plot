#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pysam
import seaborn as sns
import glob
import numpy as np
import pandas as pd
from multiprocessing import Pool

def count_mutli(bam_file):
    print 'Running %s' %bam_file.split('/')[-1].split('_')[1]
    multi_jump = 0
    aln_count = 0
    with pysam.Samfile(bam_file,'rb') as bam:
	for aln in bam:
	    if aln.is_supplementary:
		multi_jump += 1
	    aln_count += 1
    return np.true_divide(multi_jump,aln_count)

def main():
    project_path = '/scratch/02727/cdw2854/plasma_project'
    bisulfite_folder = project_path + '/bamFiles'
    bam_files =glob.glob(bisulfite_folder+'/PD_*bam')
    juming_fraction = Pool(4).map(count_mutli,bam_files)
    df = pd.DataFrame({'treatment':map(lambda x: x.split('/')[-1].split('_')[1], bam_files),
		'fractions': juming_fraction})
    print df
    sns.barplot(data = df, x = 'treatment', y ='fractions')
    plt.savefig('multi_jump_stat.pdf')

if __name__ == '__main__':
    main()
