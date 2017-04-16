from collections import defaultdict
import os
import pandas as pd
import sys

def read_fastq(file_fq):
    """
    takes a fastq file as input
    yields idSeq, sequence and score
    for each fastq entry
    learned from:
    http://codereview.stackexchange.com/questions/32897/efficient-parsing-of-fastq
    """

    cdef:
        int line_count = 0
        str line
        str idSeq, sequence, score

    while True:

        line = file_fq.readline()
        line_count += 1
        #break if we hit the end of the file
        if not line:
            break

        if line_count == 1:
            idSeq = line.strip().lstrip('@')

        elif line_count == 2:
            sequence = line.strip()

        elif line_count == 4:
            score = line.strip()
            line_count = 0
            yield idSeq,sequence, score

    yield idSeq, sequence, score

def parse_seq(sequence, seq_dict):
    cdef:
        int i
        str base

    for i, base in enumerate(sequence):
        seq_dict[base][i] += 1
    return seq_dict

def open_file(filename):
    if filename.endswith('gz') or filename.endswith('gzip'):
        return os.popen('zcat %s ' %filename)
    elif filename.endswith('fastq') or filename.endswith('fq'):
        return open(filename,'r')
    else:
        sys.exit('Unrecognized file type!')

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

    print 'Running %s' %filename
    samplename = os.path.basename(filename)
    seq_dict = defaultdict(lambda: defaultdict(int))
    with open_file(filename) as fastq:
        for fastq_count, (name, seq, qual) in enumerate(read_fastq(fastq)):
            seq_dict = parse_seq(seq, seq_dict)
            if fastq_count % 1000000 == 0:
                print 'Parsed: %i records in %s' %(fastq_count, filename)
    df = seq_dict_to_df(seq_dict) \
		.assign(sample_name = samplename) \
		.assign(position = lambda d: d.position + 1) 
    return df
