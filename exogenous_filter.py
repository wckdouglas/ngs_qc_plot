#!/usr/bin/env python

import pysam
from bwapy import BwaAligner
import argparse
import logging
import sys
import re

def getopt():
    parser = argparse.ArgumentParser(description = 'Filter read pairs with either read mapping to chrM')
    parser.add_argument('-i', '--inbam', 
                        required=True, default='-', 
                        help = 'input bam file needed to be name sorted, such that paired end reads are next to each other')
    parser.add_argument('-o','--outbam', 
                        default = '-', 
                        help = 'output bam file (defulat: - )')
    parser.add_argument('--filtered_bam', 
                        default = '', 
                        help = 'output filtered bam file (defulat: - )')
    parser.add_argument('-x', '--index',
                        required=True, 
                        help = 'Mitochondrial gnome BWA index')
    parser.add_argument('--nm',default=0.1, type=float, help='Maximum number of mismatch to tolerate')
    args = parser.parse_args()
    return args


class exogenous_mapper():
    def __init__(self, index, nm = 0.1):
        self.matched = re.compile('([0-9]+)M')
        self.clipped = re.compile('([0-9]+)S')
        self.aligner = BwaAligner(index, options = '-k 12 -B 3')
        self.NM = nm
        self.aligned=None

    def map(self, seq):
        self.aligned=None
        self.aligned = self.aligner.align_seq(seq)
        self.aligned = list(filter(self.filter_bad_cigar, self.aligned))
        return 1 if self.aligned else 0


    def filter_bad_cigar(self, aln):
        cigar = aln.cigar
        clipped_base = sum(map(int, self.clipped.findall(cigar))) or 0
        mapped_base = sum(map(int, self.matched.findall(cigar)))
        return ( (float(clipped_base)  + aln.NM)/mapped_base ) <= self.NM



def main():
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
    logger = logging.getLogger('Exogenous filter')
    args = getopt()
    pair = 0
    out_pair =0
    mapper = exogenous_mapper(args.index, nm = args.nm)
    with pysam.Samfile(args.inbam) as inbam:
        with pysam.Samfile(args.outbam, 'wb', template = inbam) as outbam:
            filter_bam = ''
            if args.filtered_bam:
                filter_bam = pysam.Samfile(args.filtered_bam, 'wb', 
                                            template = inbam)
            try:
                while True:
                    read1 = next(inbam)
                    read2 = next(inbam)
                    assert read1.query_name == read2.query_name, '%s and %s' %(read1.query_name, read2.query_name)
                    pair += 1

                    seq1 = read1.get_forward_sequence()
                    seq2 = read2.get_forward_sequence()

                    align1 = mapper.map(seq1)
                    align2 = mapper.map(seq2)

                    if align1 or align2:
                        if filter_bam:
                            filter_bam.write(read1)
                            filter_bam.write(read2)
                    else:
                        outbam.write(read1)
                        outbam.write(read2)
                        out_pair += 1
            except StopIteration:
                pass
    if filter_bam:
        filter_bam.close()
    
    logger.info('Read %i pairs, written %i pairs' %(pair, out_pair))


if __name__ == '__main__':
    main()

            
            
