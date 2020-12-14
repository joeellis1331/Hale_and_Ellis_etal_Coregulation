#!/usr/bin/python3

import argparse
import pandas
import subprocess
import requests, sys

parser = argparse.ArgumentParser(description = "conservation for Co-reg pub")
parser.add_argument("--one2one", help = "xlsx file containing the list of one-to-one orthologs between human and mouse")
parser.add_argument("--eventListMEF", help = "xlsx file of the filtered event list from rMATS")
parser.add_argument("--eventListHEK", help = "xlsx file of the filtered event list from rMATS")
parser.add_argument("--len", type = int, help = "cutoff for the length difference of exons")
args = parser.parse_args()

#read necessary files
orthoList = pandas.read_csv(args.one2one)
df_MEF = pandas.read_excel(args.eventListMEF, sheet_name = 'Raw_original', usecols = 'A:J')
df_HEK = pandas.read_excel(args.eventListHEK, sheet_name = 'Raw_original', usecols = 'A:J')


#remove genes that are not one-2-one orthologs
#use left merge to find common elements, any non-matching get a containing
#drop rows with NaN
#slice off the orthoList to clean up df a bit
df_mouse = df_MEF.merge(orthoList, how = 'left', left_on = 'GeneID', right_on = 'Mouse_ID').dropna()
#use mouse genes to look for human ones, remove extra columns from intiial comparison
df_ortho = df_HEK.merge(df_mouse, how = 'right', left_on = 'GeneID', right_on = 'Human_ID', suffixes = ('_HEK','_MEF')).dropna().iloc[:, 0:20]


#separate species back into individual datasets and drop drop_duplicates
df_ortho = df_ortho.drop_duplicates(subset=['exonStart_0base_HEK', 'exonEnd_HEK', 'exonStart_0base_MEF', 'exonEnd_MEF'])

#add empty columns to fill with genomic sequence
df_ortho.insert(loc = 10, column = 'sequence_HEK', value = '')
df_ortho.insert(loc = 21, column = 'sequence_MEF', value = '')

#convert floats to ints (needed for getting sequenes so genomic positions are 1919 not 1919.0)
df_ortho = df_ortho.astype({'exonStart_0base_HEK':'int', 'exonEnd_HEK':'int', 'upstreamES_HEK':'int',
                                    'upstreamEE_HEK':'int', 'downstreamES_HEK':'int', 'downstreamEE_HEK':'int',
                                    'exonStart_0base_MEF':'int', 'exonEnd_MEF':'int', 'upstreamES_MEF':'int',
                                    'upstreamEE_MEF':'int', 'downstreamES_MEF':'int', 'downstreamEE_MEF':'int'})



##############################################################################################
#get sequences of cassette exons
##############################################################################################

#get sequences for HEKs
for row in df_ortho.itertuples():
    #upstream exon
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/human/" + str(row.chr_HEK) + ':' + str(row.upstreamES_HEK) + '..' + str(row.upstreamEE_HEK) + ':1?'
    r_up = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r_up.ok:
      r_up.raise_for_status()
      sys.exit()

    #cassette exon
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/human/" + str(row.chr_HEK) + ':' + str(row.exonStart_0base_HEK) + '..' + str(row.exonEnd_HEK) + ':1?'
    r_cass = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r_cass.ok:
      r_cass.raise_for_status()
      sys.exit()

    #downstream exon
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/human/" + str(row.chr_HEK) + ':' + str(row.downstreamES_HEK) + '..' + str(row.downstreamEE_HEK) + ':1?'
    r_down = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r_down.ok:
      r_down.raise_for_status()
      sys.exit()

    #combine all the exonic sequences
    seq = r_up.text + r_cass.text + r_down.text
    df_ortho.at[row.Index, 'sequence_HEK'] = seq


#get sequences for MEFs
for row in df_ortho.itertuples():
    #upstream exon
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/mouse/" + str(row.chr_MEF) + ':' + str(row.upstreamES_MEF) + '..' + str(row.upstreamEE_MEF) + ':1?'
    r_up = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r_up.ok:
      r_up.raise_for_status()
      sys.exit()

    #cassette exon
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/mouse/" + str(row.chr_MEF) + ':' + str(row.exonStart_0base_MEF) + '..' + str(row.exonEnd_MEF) + ':1?'
    r_cass = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r_cass.ok:
      r_cass.raise_for_status()
      sys.exit()

    #downstream exon
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/mouse/" + str(row.chr_MEF) + ':' + str(row.downstreamES_MEF) + '..' + str(row.downstreamEE_MEF) + ':1?'
    r_down = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

    if not r_down.ok:
      r_down.raise_for_status()
      sys.exit()

    #combine all the exonic sequences
    seq = r_up.text + r_cass.text + r_down.text
    df_ortho.at[row.Index, 'sequence_MEF'] = seq



##############################################################################################################
#compare length and remove any sequences that are not of similar length
#############################################################################################################
#df_ortho = df_ortho[abs(df_ortho['sequence_HEK'].str.len() - df_ortho['sequence_MEF'].str.len()) <= args.len]


###############################################################################################################
#write to a fasta record
###############################################################################################################
h = open('ortho_fasta_HEK.txt', 'w+')
m = open('ortho_fasta_MEF.txt', 'w+')

for row in df_ortho.itertuples():
    h.write('>' + row.geneSymbol_HEK + '_' + row.chr_HEK+ ':' + str(row.upstreamES_HEK) + '..' + str(row.downstreamEE_HEK) + '\n' + row.sequence_HEK + '\n')
    m.write('>' + row.geneSymbol_MEF + '_' + row.chr_MEF+ ':' + str(row.upstreamES_MEF) + '..' + str(row.downstreamEE_MEF) + '\n' + row.sequence_MEF + '\n')

h.close
m.close
