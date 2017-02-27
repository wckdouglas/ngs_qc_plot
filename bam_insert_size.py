
#!/bin/env python

import sys
if len(sys.argv) != 3:
    sys.exit('[usage] python %s <bam_file> <out_file>' %(sys.argv[0]))

import numpy as np
import pysam
import glob
import pandas as pd
import os

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
    df = pd.DataFrame({'isize': fragSize, 'counts':count}) \
        .assign(percentage = lambda d: np.true_divide(d['counts'], d['counts'].sum())) \
        .assign(samplename = samplename)
    return df


def main():
    bam_file = sys.argv[1]
    out_file = sys.argv[2]
    if not os.path.isfile(bam_file):
        sys.exit('[usage] python %s <bam_file> <out_file>' %(sys.argv[0])+\
                '[ERROR] %s is not a directory' %datapath)

    df = parseBam(bam_file)
    df.to_csv(out_file,sep='\t',index=False)
    return 0


if __name__ == '__main__':
    main()
