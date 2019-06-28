#!/usr/bin/env python

import pysam
import re
from sequencing_tools.bam_tools import cigar_to_str
from sequencing_tools.io_tools import xopen

soft_clipped = re.compile('S+')
clipped5 = 'clipped5.fa'
clipped3 = 'clipped5.fa'
clipped5_count = 0
clipped3_count = 0
with pysam.Samfile('-', 'rb') as inbam, \
        xopen(clipped5,'w') as clip5, \
        xopen(clipped3, 'w') as clip3:
    for aln in inbam:
        if not aln.is_unmapped and 'S' in aln.cigarstring:
            cigar = cigar_to_str(aln.cigarstring)
            for m in soft_clipped.finditer(cigar):
                if m.start() == 0:
                    print('>{read_id}\n{seq}'.format(read_id=aln.query_name,
                                                     seq = aln.query_sequence[:m.end()]), 
                          file = clip5)
                    clipped5_count += 1
                else:
                    print('>{read_id}\n{seq}'.format(read_id=aln.query_name,
                                                     seq = aln.query_sequence[m.start():m.end()]),
                          file = clip3)
                    clipped3_count += 1
print("Extracted %i 5' softclip and %i 3' softclip" %(clipped5_count, clipped3_count))
