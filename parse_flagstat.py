#!/usr/bin/env python

import glob
import sys
import os
import pandas as pd

def get_number(line):
    return [int(line.split(' ')[0])]

def get_line(lines, keyword):
    return_line = ''
    for line in lines:
        if keyword in line:
            return_line += line 
    return return_line


def parse_stat(stat_file):
    samplename = os.path.basename(stat_file)
    info = open(stat_file,'r').readlines()
    stat_df = pd.DataFrame()
    stat_df['duplicates'] = get_number(get_line(info, 'duplicates'))
    stat_df['trimmed reads'] = get_number(get_line(info, 'read1'))
    stat_df['mapped'] = get_number(get_line(info, 'mapped'))
    stat_df['supplementary'] = get_number(get_line(info, 'supplementary'))
    stat_df['proper pair'] = get_number(get_line(info, 'properly paired'))
    stat_df['samplename'] = samplename
    return stat_df
    

def main():
    if len(sys.argv) != 2:
        sys.exit('[usage] python %s <flagstat_files_path>\n' %sys.argv[0])
    datapath = sys.argv[1]
    tablename = datapath + '/combined_stats.tsv' 
    stat_files = glob.glob(datapath + '/*.stats')
    df = map(parse_stat,stat_files)
    df = pd.concat(df)
    df.to_csv(tablename, sep='\t',index=False)
    print 'Written: %s' %tablename

if __name__ == '__main__':
    main()
