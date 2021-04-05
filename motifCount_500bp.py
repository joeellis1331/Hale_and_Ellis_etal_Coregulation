#!/usr/bin/python3

import argparse
import numpy
import pandas
from Bio import SeqIO
import re
from random import shuffle
from collections import Counter


#args
parser = argparse.ArgumentParser(description="motif stuff for CoReg paper")
parser.add_argument("--fasta", help="<Cell>eventsSCORE.fa")
parser.add_argument("--SCORE", help="SCORE spreadsheet from MBNLsig_common_dd<cell>")
args = parser.parse_args()

###########set the stage and import files###############################

#block sizes to differentiate exon/intron sections based on full seq length
dfInc = pandas.read_excel(args.SCORE, sheet_name = "Reorg + SCORE Inclu", usecols = "B, G, J, K, T")
dfExc = pandas.read_excel(args.SCORE, sheet_name = "Reorg + SCORE Exclu", usecols = "B, G, J, K, T")
dfSCORE = pandas.concat([dfInc, dfExc])
dfSCORE['upstreamES'] = dfSCORE['upstreamES'] + 1 #NOTE: For whatever reason the fasta grabbed from UCSC and my SCORE MBNLsig common files do not match coordinates. To circumvent this I added one nucleotide to the upstreamES
dfSCORE['indexer'] = dfSCORE['geneSymbol'] + ':' + dfSCORE['upstreamES'].astype(str) + '-' + dfSCORE['downstreamEE'].astype(str) #creates matching indexer
dfSCORE = dfSCORE[['indexer', 'IncLevelDifference_d0r0_d1000r0', 'SCORE_0.1']] #rearrange columns


#empty dataframe to populate
df500 = pandas.DataFrame(columns = ['indexer', 'geneSymbol', 'upIntronLength', 'downIntronLength',
                        'upIntron nt ratio (A, G, C, T)', 'downIntron nt ratio (A, G, C, T)',
                        'YGCY_up', 'YGCY_down',
                        'GCAUG_up', 'GCAUG_down',
                        'UGCAUG_up', 'UGCAUG_down', 'AGCAUG_up', 'AGCAUG_down', 'GGCAUG_up', 'GGCAUG_down', 'CGCAUG_up', 'CGCAUG_down',
                        'RGCAUG_up', 'RGCAUG_down', 'RGCGUG_up', 'RGCGUG_down', 'RGCCUG_up', 'RGCCUG_down', 'RGCUUG_up', 'RGCUUG_down',
                        'YGCAUG_up', 'YGCAUG_down', 'YGCGUG_up', 'YGCGUG_down', 'YGCCUG_up', 'YGCCUG_down', 'YGCUUG_up', 'YGCUUG_down'])


##################populate the spreadsheet and find motifs and etc#######################
#function to parse a string into kmers
def get_k_mer(string, k):
   length = len(string)
   return [string[i: i+ k] for i in range(length-k+1)]

############################### meat n potata's ######################################################
#read in fasta record from getFasta outputfile using SeqIO
with open(args.fasta) as handle:
    for record in SeqIO.FastaIO.FastaIterator(handle): #note: this utilizes a SeqRecord class
        #use str() to convert the SeqRecord class obj into useable str
        name = (str(record.id)).split('_')[4] #grabs the gene symbol
        description = (str(record.description)).split(' ') #grabs the entire fasta header for coordinate indexing
        chrCoord = description[1].split(':')[1] #grabs the start and end coordinates
        indexer = name + ':' + chrCoord #combines gene symbol with start and end coordinates to be able to index

        ############################### grab sequence information and split, scramble, limit ############################################
        #splits sequence into exon and intron parts
        seq = re.sub(r'((?<=[a-z])[A-Z]|(?<!\A)[A-Z](?=[a-z]))', r' \1', str(record.seq))  #split based on upper and lower case str(record.seq)
        seq = seq.split()
        #apply separate chunks to appropriate intron/exon labels
        upExon = seq[0]
        upIntron = seq[1]
        upIntron = upIntron.lower()
        upIntronLength = len(upIntron)
        cassExon = seq[2]
        downIntron = seq[3]
        downIntron = downIntron.lower()
        downIntronLength = len(downIntron)
        downExon = seq[4]


        #grab a window of intronic space
        upIntron500 = upIntron[-501:-1]
        upIntron500Length = len(upIntron500)
        downIntron500 = downIntron[:500]
        downIntron500Length = len(downIntron500)

        #sequence composition of intronic space
        upA, upG, upC, upT = (upIntron500.count('a')/500) , (upIntron500.count('g')/500), (upIntron500.count('c')/500), (upIntron500.count('t')/500)
        upNT = [upA, upG, upC, upT]
        downA, downG, downC, downT = (downIntron.count('a')/500), (downIntron.count('g')/500), (downIntron.count('c')/500), (downIntron.count('t')/500)
        downNT = [downA, downG, downC, downT]

        #using Counter function to find the number of each kmer present in the intron/exon
        motifs_upIntron_4 = Counter(get_k_mer(upIntron500, 4))
        motifs_downIntron_4 = Counter(get_k_mer(downIntron500, 4))
        motifs_upIntron_5 = Counter(get_k_mer(upIntron500, 5))
        motifs_downIntron_5 = Counter(get_k_mer(downIntron500, 5))
        motifs_upIntron_6 = Counter(get_k_mer(upIntron500, 6))
        motifs_downIntron_6 = Counter(get_k_mer(downIntron500, 6))
        motifs_upIntron_7 = Counter(get_k_mer(upIntron500, 7))
        motifs_downIntron_7 = Counter(get_k_mer(downIntron500, 7))

        #count the number of YGCYs and UGCAUGs
        #MBNL motifs
        YGCY_up = motifs_upIntron_4['cgcc'] + motifs_upIntron_4['tgcc'] + motifs_upIntron_4['cgct'] + motifs_upIntron_4['tgct']
        YGCY_down = motifs_downIntron_4['cgcc'] + motifs_downIntron_4['tgcc'] + motifs_downIntron_4['cgct'] + motifs_downIntron_4['tgct']

        # GCTTGCT_up, CGCTTGC_up, GCTGCTT_up, TGCTTGC_up, GCTTCGC_up  = motifs_upIntron_7['gcttgct'], motifs_upIntron_7['cgcttgc'], motifs_upIntron_7['gctgctt'], motifs_upIntron_7['tgcttgc'], motifs_upIntron_7['gcttcgc']
        # GCTTGCT_down, CGCTTGC_down, GCTGCTT_down, TGCTTGC_down, GCTTCGC_down  = motifs_downIntron_7['gcttgct'], motifs_downIntron_7['cgcttgc'], motifs_downIntron_7['gctgctt'], motifs_downIntron_7['tgcttgc'], motifs_downIntron_7['gcttcgc']

        #RBFOX motifs
        GCAUG_up = motifs_upIntron_5['gcatg']
        GCAUG_down = motifs_downIntron_5['gcatg']
        UGCAUG_up = motifs_upIntron_6['tgcatg']
        UGCAUG_down = motifs_downIntron_6['tgcatg']
        # GCAUGC_up, GCAUGU_up, AGCAUG_up, GCAUGA_up = motifs_upIntron_6['gcatgc'], motifs_upIntron_6['gcatgt'], motifs_upIntron_6['agcatg'], motifs_upIntron_6['gcatga']
        # GCAUGC_down, GCAUGU_down, AGCAUG_down, GCAUGA_down = motifs_downIntron_6['gcatgc'], motifs_downIntron_6['gcatgt'], motifs_downIntron_6['agcatg'], motifs_downIntron_6['gcatga']

        #NGCAUG motifs (already have UGCAUG)
        AGCAUG_up, AGCAUG_down = motifs_upIntron_6['agcatg'], motifs_downIntron_6['agcatg']
        GGCAUG_up, GGCAUG_down = motifs_upIntron_6['ggcatg'], motifs_downIntron_6['ggcatg']
        CGCAUG_up, CGCAUG_down = motifs_upIntron_6['cgcatg'], motifs_downIntron_6['cgcatg']

        #specific SCORE motifs
        #RGCNUG
        RGCAUG_up = motifs_upIntron_6['agcatg'] + motifs_upIntron_6['ggcatg']
        RGCAUG_down = motifs_downIntron_6['agcatg'] + motifs_downIntron_6['ggcatg']
        RGCGUG_up = motifs_upIntron_6['agcgtg'] + motifs_upIntron_6['ggcgtg']
        RGCGUG_down = motifs_downIntron_6['agcgtg'] + motifs_downIntron_6['ggcgtg']
        RGCCUG_up = motifs_upIntron_6['agcctg'] + motifs_upIntron_6['ggcctg']
        RGCCUG_down = motifs_downIntron_6['agcctg'] + motifs_downIntron_6['ggcctg']
        RGCUUG_up = motifs_upIntron_6['agcttg'] + motifs_upIntron_6['ggcttg']
        RGCUUG_down = motifs_downIntron_6['agcttg'] + motifs_downIntron_6['ggcttg']

        YGCAUG_up = motifs_upIntron_6['cgcatg'] + motifs_upIntron_6['tgcatg']
        YGCAUG_down = motifs_downIntron_6['cgcatg'] + motifs_downIntron_6['tgcatg']
        YGCGUG_up = motifs_upIntron_6['cgcgtg'] + motifs_upIntron_6['tgcgtg']
        YGCGUG_down = motifs_downIntron_6['cgcgtg'] + motifs_downIntron_6['tgcgtg']
        YGCCUG_up = motifs_upIntron_6['cgcctg'] + motifs_upIntron_6['tgcctg']
        YGCCUG_down = motifs_downIntron_6['cgcctg'] + motifs_downIntron_6['tgcctg']
        YGCUUG_up = motifs_upIntron_6['cgcttg'] + motifs_upIntron_6['tgcttg']
        YGCUUG_down = motifs_downIntron_6['cgcttg'] + motifs_downIntron_6['tgcttg']


        df500 = df500.append({'indexer':indexer, 'geneSymbol':name, 'upIntronLength':upIntron500Length, 'downIntronLength':downIntron500Length,
                                'upIntron nt ratio (A, G, C, T)':upNT, 'downIntron nt ratio (A, G, C, T)':downNT,
                                'YGCY_up': YGCY_up, 'YGCY_down':YGCY_down,
                                'GCAUG_up':GCAUG_up, 'GCAUG_down':GCAUG_down,
                                'UGCAUG_up':UGCAUG_up, 'UGCAUG_down':UGCAUG_down, 'AGCAUG_up':AGCAUG_up, 'AGCAUG_down':AGCAUG_down, 'GGCAUG_up':GGCAUG_up, 'GGCAUG_down':GGCAUG_down, 'CGCAUG_up':CGCAUG_up, 'CGCAUG_down':CGCAUG_down,
                                'RGCAUG_up':RGCAUG_up, 'RGCAUG_down':RGCAUG_down,'RGCGUG_up':RGCGUG_up, 'RGCGUG_down':RGCGUG_down, 'RGCCUG_up':RGCCUG_up, 'RGCCUG_down':RGCCUG_down, 'RGCUUG_up':RGCUUG_up, 'RGCUUG_down':RGCUUG_down,
                                'YGCAUG_up':YGCAUG_up, 'YGCAUG_down':YGCAUG_down,'YGCGUG_up':YGCGUG_up, 'YGCGUG_down':YGCGUG_down, 'YGCCUG_up':YGCCUG_up, 'YGCCUG_down':YGCCUG_down, 'YGCUUG_up':YGCUUG_up, 'YGCUUG_down':YGCUUG_down}, ignore_index = True)





df500 = df500.merge(dfSCORE, how = 'outer', on = 'indexer')
df500 = df500.drop_duplicates(subset = ['geneSymbol', 'IncLevelDifference_d0r0_d1000r0'], keep = 'first')

df500.to_excel('cisCountsHEK_500bp.xlsx', index = False)





#############################
# df500 = df500.append({'indexer':indexer, 'geneSymbol':name, 'upIntronLength':upIntron500Length, 'downIntronLength':downIntron500Length,
#                         'upIntron nt ratio (A, G, C, T)':upNT, 'downIntron nt ratio (A, G, C, T)':downNT,
#                         'YGCY_up': YGCY_up, 'YGCY_down':YGCY_down, 'GCTTGCT_up':GCTTGCT_up, 'GCTTGCT_down':GCTTGCT_down, 'CGCTTGC_up':CGCTTGC_up, 'CGCTTGC_down':CGCTTGC_down, 'GCTGCTT_up':GCTGCTT_up, 'GCTGCTT_down':GCTGCTT_down, 'TGCTTGC_up':TGCTTGC_up, 'TGCTTGC_down':TGCTTGC_down, 'GCTTCGC_up':GCTTCGC_up, 'GCTTCGC_down':GCTTCGC_down,
#                         'GCAUG_up':GCAUG_up, 'GCAUG_down':GCAUG_down, 'UGCAUG_up':UGCAUG_up, 'UGCAUG_down':UGCAUG_down, 'GCAUGC_up':GCAUGC_up, 'GCAUGC_down':GCAUGC_down, 'GCAUGU_up':GCAUGU_up, 'GCAUGU_down':GCAUGU_down, 'AGCAUG_up':AGCAUG_up, 'AGCAUG_down':AGCAUG_down, 'GCAUGA_up':GCAUGA_up, 'GCAUGA_down':GCAUGA_down,
#                         'RGCAUG_up':RGCAUG_up, 'RGCAUG_down':RGCAUG_down, 'YGCAUG_up':YGCAUG_up, 'YGCAUG_down':YGCAUG_down}, ignore_index = True)
#
