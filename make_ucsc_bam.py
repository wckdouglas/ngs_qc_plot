#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import stat
import shutil

def errorRemoveReadonly(func, path, exc):
    excvalue = exc[1]
    if func in (os.rmdir, os.remove) and excvalue.errno == errno.EACCES:
        # change the file to be readable,writable,executable: 0777
        os.chmod(path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)  
        # retry
        func(path)
    else:
        # raiseenter code here

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
shutil.rmtree(tmp_dir, ignore_errors=False, onerror=errorRemoveReadonly) 
