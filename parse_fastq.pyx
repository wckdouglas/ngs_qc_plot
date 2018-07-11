from collections import defaultdict
import os
import pandas as pd
import sys
from sequencing_tools.io_tools import xopen
from sequencing_tools.fastq_tools import readfq

def parse_seq(sequence, seq_dict):
    cdef:
        int i
        str base

    for i, base in enumerate(sequence):
        seq_dict[base][i] += 1
    return seq_dict


def seq_dict_to_df(seq_dict):
    d = []
    for base, pos_counts in seq_dict.iteritems():
        d.append(
                pd.DataFrame({'position':pos_counts.keys(),
                            'count':pos_counts.values()}) \
                    .assign(base = base)
                )
    return pd.concat(d, axis=0)


def read_file(filename):
    cdef:
        int fastq_count = 0
        str name, seq, qual

    print('Running %s' %filename)
    samplename = os.path.basename(filename)
    seq_dict = defaultdict(lambda: defaultdict(int))
    with xopen(filename, 'r') as fastq:
        for fastq_count, (name, seq, qual) in enumerate(readfq(fastq)):
            seq_dict = parse_seq(seq, seq_dict)
            if fastq_count % 1000000 == 0:
                print('Parsed: %i records in %s' %(fastq_count, filename))
    df = seq_dict_to_df(seq_dict) \
        .assign(sample_name = samplename) \
        .assign(position = lambda d: d.position + 1) 
    return df
