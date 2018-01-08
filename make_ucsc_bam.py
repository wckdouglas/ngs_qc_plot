#!/usr/bin/env python

from __future__ import print_function
import os
import sys

if len(sys.argv) < 3:
    sys.exit('[usage] python %s <inbam> <outbam>' %(sys.argv[0]))

in_bam = sys.argv[1]
out_bam = sys.argv[2]
threads = int(sys.argv[3])

tmp_dir = out_bam.replace('.bam','')
if not os.path.isdir(tmp_dir):
    os.mkdir(tmp_dir)

command = 'samtools view -h@ %i %s ' %(threads, in_bam)+\
        "| grep -v 'gi|23898|emb|X12811.1|' "+\
        "| grep -v 'gi|555853|gb|U13369.1|HSU13369' " +\
        "| samtools view -b@ %i " %(threads)+\
        "| sambamba sort --show-progress --nthreads %i "%(threads)+\
        " -o %s --tmpdir=%s /dev/stdin" %(out_bam, tmp_dir) 
print('Running: ', command, file=sys.stderr)
os.system(command)
os.removedirs(tmp_dir)
