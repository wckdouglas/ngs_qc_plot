#!/usr/bin/env python

from glob import glob 
import re
import numpy as np
import os
import argparse


def getopt():
    parser = argparse.ArgumentParser(description = 'Merging files with similar name [01-$samplename.[bam|bed]]')
    parser.add_argument('-p', '--path', required=True, help = 'file path')
    parser.add_argument('-s', '--suffix', required=True, choices = ['fastq.gz','fastq','bed','bam'], help = 'file suffix')
    parser.add_argument('-d','--dry', action='store_true', help = 'print out commands to be run')
    parser.add_argument('-r','--rmdup', action='store_true', help ='Only work for bed file')
    args = parser.parse_args()
    return args

def get_group(y):
    return map(lambda x: re.sub('^[0-9]+-','',os.path.basename(x)), y)

def make_dir(directory):
    if not os.path.isdir(directory):
        os.mkdir(directory)

def main():
    args = getopt()
    path = args.path
    outpath = path + '/merged' if not args.rmdup else path + '/merged_rmdup'
    file_type = args.suffix
    make_dir(outpath)


    map_files = np.array(glob(path + '/*' + file_type))
    group_name = np.array(get_group(map_files))

    group_file = np.unique(group_name)

    for group in group_file:
        files = map_files[group_name==group]
        joined_files = ' '.join(files)
        if file_type == 'bam':
            command = 'samtools cat %s > %s/%s ' %(joined_files, outpath, group)
        elif not args.rmdup:
            command = 'cat %s > %s/%s ' %(joined_files, outpath, group)
        elif args.rmdup:
            if file_type == 'bed':
                temp_dir = ('%s/%s' %(outpath, group)).replace('.bed','')
                make_dir(temp_dir)            
                command = 'cat %s ' %(joined_files) + \
                    '| sort -k1,1 -k2,2n -k3,3n -k6,6 '+\
                        '--temporary-directory=%s -u ' %(temp_dir)  +\
                    '> %s/%s ' %(outpath, group)
            else:
                sys.exit('Cannot run rmdup with files other than BED')
        if not args.dry:
            os.system(command)
        print command

if __name__ == '__main__':
    main()
