#!/usr/bin/env python

import sys
if len(sys.argv) < 2:
    sys.exit('[usage] python %s <fastqc_path1> [fastqc_path2] [fastqc_path3] ... ')

from bs4 import BeautifulSoup
import pandas as pd
import os
import glob

def read_fastqc(html_file):
    samplename = os.path.basename(html_file)
    with open(html_file) as f:
        html = BeautifulSoup(f, "lxml")
        tab = html.find_all('table')
        df = pd.read_html(str(tab[0]))[0] \
            .assign(samplename = samplename)

    return df


def main():
    fastqc_paths = sys.argv[1:]
    dfs = []

    for fp in fastqc_paths:
        htmls = glob.glob(fp + '/*fastqc.html')
        for html in htmls:
            df = read_fastqc(html)
            dfs.append(df)

    df = pd.concat(dfs) \
        .pipe(pd.pivot_table, index = 'samplename',
              values = 'Value', columns = 'Measure',
              aggfunc=lambda x: ','.join(x)) 
    df.to_csv('fastqc_profile.tsv', sep='\t', index=False)



if __name__ == '__main__':
    main()
