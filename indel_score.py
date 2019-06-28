#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
import sys
from multiprocessing import Pool

def generate_index_array(sequence):
    current_base = 'X'
    base_count = 0

    repeat_index = []
    for i, base in enumerate(sequence):
        if base == current_base:
            base_count += 1
        else:
            current_base = base
            base_count = 0
        repeat_index.append(base_count)
        if i % 1000000 == 0:
            print('Parsed %i position' %i, file = sys.stderr)
    return repeat_index

def main():
    if len(sys.argv) != 2:
        sys.exit('[usage]: python %s <fasta file>' %sys.argv[0])


    for i, record in enumerate(SeqIO.parse(sys.argv[1],'fasta')):
        sequence = record.seq
        if i == 0:
            mode = 'w'
            header = True

        else:
            mode = 'a'
            header=False

        positive_index, negative_index = Pool(2).map(generate_index_array, [sequence, sequence.reverse_complement()])
        df = pd.DataFrame({'start':range(len(sequence)),
                        'positive_index':positive_index,
                        'negative_index':negative_index[::-1]})\
            .assign(fwd_base = list(str(sequence))) \
            .assign(contig = record.id)
        df.to_csv(sys.stdout,
                  sep='\t',index=False,
                  header=header,
                  mode = mode)
        print('Finished %s' %record.id, file = sys.stderr)


if __name__ == '__main__':
    main()
