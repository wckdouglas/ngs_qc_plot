#!usr/bin/env python


import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Converting csv files to excels, split column into sheets')
parser.add_argument('--csv', required=True,  help='input csv file')
parser.add_argument('--excel', required=True, help ='output excel file')
parser.add_argument('--sheet', default=None, help ='Sheet name, must be a column name in the csv file')
args = parser.parse_args()

df = pd.read_csv(args.csv)
writer = pd.ExcelWriter(args.excel)
if not args.sheet:
    df.to_excel(writer,'Sheet1')

else:
    assert(args.sheet in df.columns):
        for col, col_data in df.groupby(args.sheet):
            col_data.to_excel(writer, col)
            print('Written %s' %col)
print('Written %s' %args.excel)




