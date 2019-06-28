#!/usr/bin/env python

import sys
from operator import itemgetter
from collections import defaultdict

#gtf = '/stor/work/Lambowitz/ref/benchmarking/GRCH38_genome/genes.gtf'
if len(sys.argv) != 2:
    sys.exit('[usage] python %s <gtf file>' %sys.argv[0])
gtf = sys.argv[1]


def parse_extra_fields(extra_field):
    info_fields = extra_field.split(';')
    info_dict = defaultdict(str)
    for info in info_fields[:-1]:
        row_fields = info.strip().split('=') if '=' in info else info.strip().split(' ')
        info_dict[row_fields[0]] = row_fields[1].strip('"')

    return info_dict



if __name__ == '__main__':
    for line in open(gtf,'r'):
        if not line.startswith('#'):
            fields = line.strip().split('\t')
            if fields[2] == "gene":
                chrom, start, end, \
                strand, extra_fields = itemgetter(0,3,4,6,8)(fields)
                info_dict = parse_extra_fields(extra_fields)


                line = '{chrom}\t{start}\t{end}\t{gene_name}\t'\
                        '0\t{strand}\t{gene_type}\t{gene_id}' \
                        .format(chrom = chrom,
                                start = start, 
                                end = end,
                                strand = strand,
                                gene_name = info_dict['gene_name'] or info_dict['Name'] or info_dict['gene_id'],
                                gene_id = info_dict['gene_id'],
                                gene_type = info_dict['gene_biotype'] or info_dict['biotype'])
                print(line, file = sys.stdout)
