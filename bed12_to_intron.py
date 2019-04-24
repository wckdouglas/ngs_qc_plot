#!/usr/bin/env python

import fileinput
from operator import itemgetter
import numpy as np
from sequencing_tools.gene_tools import Bed12Record

'''
usage: cat bed12 | python bed12_to_intron.py  > intron.bed
'''


def main():
    for line in fileinput.input():
        gr = Bed12Record(line)
        for intron in gr.get_introns():
            print(intron)



if __name__ == '__main__':
    main()
        
