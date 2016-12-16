#!/usr/bin/env python

from Bio import SeqIO
import sys
import argparse



def parse_opt():
    parser = argparse.ArgumentParser(description='Convert genome fasta file to other format.')
    parser.add_argument('-f','--fasta', required=True,
                        help = 'Input genome fasta file')
    parser.add_argument('-o','--out_fmt', choices = ['bed','genome','header'], required=True,
                        help = 'Output format')
    return parser.parse_args()


def make_sam_header(fasta_file):
    for record in SeqIO.parse(fasta_file,'fasta'):
        print '@SQ\tSN:%s\tLN:%i' %(record.id, len(record.seq))

def make_genome(fasta_file):
    for record in SeqIO.parse(fasta_file,'fasta'):
        print '%s\t%i' %(record.id, len(record.seq))
        
def make_bed(fasta_file):
    for record in SeqIO.parse(fasta_file,'fasta'):
        print '%s\t1\t%i\t%s\t0\t+' %(record.id, len(record.seq),
                                    record.id)
        
def main():
    args = parse_opt()
    if args.out_fmt == 'genome':
        make_genome(args.fasta)
    elif args.out_fmt == 'bed':
        make_bed(args.fasta)
    elif args.out_fmt == 'header':
        make_sam_header(args.fasta)

if __name__ == '__main__':
    main()