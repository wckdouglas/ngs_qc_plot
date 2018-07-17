#!/usr/bin/env python3

import pysam
import sys 
import re
from sequencing_tools.fastq_tools import reverse_complement


if 3 > len(sys.argv) or len(sys.argv) > 4:
    sys.exit('[usage] python %s <fa_file> <chr:start-end> [strand]' %sys.argv[0])

fa = sys.argv[1]
region = sys.argv[2]
strand = '+'
try:
    strand = sys.argv[3]
except IndexError:
    pass
region = re.sub(',','', region)
matches = re.search('(chr[0-9XY]+):([0-9]+)-([0-9]+)', region)
chrom, start, end = matches.groups()

fa = pysam.Fastafile(fa)
seq = fa.fetch(chrom, int(start), int(end))


seq = reverse_complement(seq) if strand == "-" or strand=="reverse" else seq

print('>' + region + ';' +strand)
print(seq)
