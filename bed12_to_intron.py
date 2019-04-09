#!/usr/bin/env python

import fileinput
from operator import itemgetter


'''
usage: cat bed12 | python bed12_to_intron.py  > intron.bed
'''


class GeneRecord():
    def __init__(self, line):
        fields = line.strip().split('\t')
        chrom, start, end, tid, exon_count, \
            strand, exon_starts, exon_ends = itemgetter(0,1,2,3,4,5,8,9)(fields)
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.exon_count = int(exon_count)
        self.tid = tid
        self.exon_starts = list(map(int, exon_starts.strip(',').split(',')))
        self.exon_ends = list(map(int, exon_ends.strip(',').split(',')))
        assert(len(self.exon_starts)==len(self.exon_ends))

    def get_introns(self):
        for i, (next_s, end) in enumerate(zip(self.exon_starts[1:], self.exon_ends)):
            intron_count = self.exon_count - i if self.strand == '-' else i + 1
            intron_line = '{chrom}\t{start}\t{end}\t{tid}\t{intron_count}\t{strand}'\
                    .format(chrom = self.chrom,
                            start = end,
                            end = next_s,
                            tid = self.tid,
                            intron_count = intron_count,
                            strand = self.strand)
            yield intron_line



def main():
    for line in fileinput.input():
        gr = GeneRecord(line)
        for intron in gr.get_introns():
            print(intron)



if __name__ == '__main__':
    main()
        
