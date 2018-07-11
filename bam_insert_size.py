#!/bin/env python


from __future__ import print_function
import sys
if 3 > len(sys.argv) or len(sys.argv) > 4:
    sys.exit('[usage] python %s <bam_file> <out_file> <mq>' %(sys.argv[0]))

import numpy as np
import pysam
import glob
import pandas as pd
import os

def filter_align(aln, mq):
    if aln.is_proper_pair and not aln.is_secondary \
            and aln.is_read1 and not aln.is_supplementary \
            and not aln.is_duplicate and aln.mapping_quality >= mq:
        return abs(aln.isize)
    else:
        return 0

def parseBam(bamFile, mq):
    samplename = os.path.basename(bamFile).split('.')[0]
    print('Extracting %s' %samplename, file = sys.stderr)
    with pysam.Samfile(bamFile, 'rb') as bam:
        fragSize = np.abs([filter_align(aln, mq) for aln in bam])
    fragSize = fragSize[fragSize>20]
    fragSize, count = np.unique(fragSize, return_counts=True)
    df = pd.DataFrame({'isize': fragSize, 'counts':count}) \
        .assign(percentage = lambda d: np.true_divide(d['counts'], d['counts'].sum())) \
        .assign(samplename = samplename)
    return df


def main():

    bam_file = sys.argv[1]
    out_file = sys.argv[2]
    mq = int(sys.argv[3]) or 0

    if not os.path.isfile(bam_file):
        sys.exit('[usage] python %s <bam_file> <out_file> <mq>' %(sys.argv[0])+\
                '[ERROR] %s is not a directory' %datapath)

    print('Using MAPQ >= %i' %(mq))
    df = parseBam(bam_file, mq)
    df.to_csv(out_file,sep='\t',index=False)
    return 0


if __name__ == '__main__':
    main()
