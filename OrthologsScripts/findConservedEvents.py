#!/usr/bin/python3

import argparse
import pandas
import re

parser = argparse.ArgumentParser(description = "conservation for Co-reg pub")
parser.add_argument("--conEvents", help = "xlsx file containing the list of one-to-one orthologs between human and mouse")
parser.add_argument("--MEF", help = "xlsx file of the filtered event list from rMATS")
parser.add_argument("--HEK", help = "xlsx file of the filtered event list from rMATS")
args = parser.parse_args()

con_MEF = pandas.read_csv(args.conEvents)
con_HEK = con_MEF
df_MEF = pandas.read_excel(args.MEF, sheet_name = 'Reorganized')
df_HEK = pandas.read_excel(args.HEK, sheet_name = 'Reorganized')


#add columns to the ConEvents df and grab appropriate string slices. NOTE: need surrounding parantheses to mark capture groups
#index locations for MEFs
con_MEF['geneSymbol'] = con_MEF['qseqid'].str.extract(r'([^_]*)')
con_MEF['upstreamES'] = con_MEF['qseqid'].str.extract(r'\:(.*?)\.')
con_MEF['downstreamEE'] = con_MEF['qseqid'].str.extract(r'([^.]+$)')
con_MEF = con_MEF.astype({'geneSymbol':'str', 'upstreamES':'int', 'downstreamEE':'int'})
df_MEF = df_MEF.merge(con_MEF, how = 'left', on = ['geneSymbol', 'upstreamES', 'downstreamEE'])
df_MEF.to_csv("MEFout.csv")


con_HEK['geneSymbol'] = con_HEK['sseqid'].str.extract(r'([^_]*)')
con_HEK['upstreamES'] = con_HEK['sseqid'].str.extract(r'\:(.*?)\.')
con_HEK['downstreamEE'] = con_HEK['sseqid'].str.extract(r'([^.]+$)')
con_HEK = con_HEK.astype({'geneSymbol':'str', 'upstreamES':'int', 'downstreamEE':'int'})
df_HEK = df_HEK.merge(con_HEK, how = 'left', on = ['geneSymbol', 'upstreamES', 'downstreamEE'])
df_HEK.to_csv('HEKout.csv')
