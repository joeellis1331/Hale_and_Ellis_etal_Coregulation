#!/usr/bin/python3

import sys
import numpy #pandas is a dependeny of numpy
import pandas

#NOT .read_csv that uses a standard '' separator. read_table uses '\t'
df_MBNL = pandas.read_table(sys.argv[1]) #reads d0r0_d1000r0 SE.JCEC
df_FOX = pandas.read_table(sys.argv[2]) #reads d0r0_d0r1 SE.JCEC
df_MF = pandas.read_table(sys.argv[3]) #reads d0r0_d1000r1 SE.JCEC
print('df_MBNL\n',df_MBNL, '\n', 'df_FOX\n', df_FOX, '\n', 'df_MF\n', df_MF, '\n')

#rename common column names to ensure there is no mixing of data (delta PSI etc)
df_MBNL = df_MBNL.rename({"IJC_SAMPLE_1":"IJC_d0r0","SJC_SAMPLE_1":"SJC_d0r0",
                            "IJC_SAMPLE_2":"IJC_d1000r0","SJC_SAMPLE_2":"SJC_d1000r0",
                            "PValue":"PValue_d0r0_d1000r0", "FDR":"FDR_d0r0_d1000r0","IncLevel1":"IncLevel_d0r0",
                            "IncLevel2":"IncLevel_d1000r0","IncLevelDifference":"IncLevelDifference_d0r0_d1000r0"}, axis='columns')

df_FOX = df_FOX.rename({"IJC_SAMPLE_1":"IJC_d0r0","SJC_SAMPLE_1":"SJC_d0r0",
                            "IJC_SAMPLE_2":"IJC_d0r1","SJC_SAMPLE_2":"SJC_d0r1",
                            "PValue":"PValue_d0r0_d0r1", "FDR":"FDR_d0r0_d0r1","IncLevel1":"IncLevel_d0r0",
                            "IncLevel2":"IncLevel_d0r1","IncLevelDifference":"IncLevelDifference_d0r0_d0r1"}, axis='columns')

df_MF = df_MF.rename({"IJC_SAMPLE_1":"IJC_d0r0","SJC_SAMPLE_1":"SJC_d0r0",
                            "IJC_SAMPLE_2":"IJC_d1000r1","SJC_SAMPLE_2":"SJC_d1000r1",
                            "PValue":"PValue_d0r0_d1000r1","FDR":"FDR_d0r0_d1000r1","IncLevel1":"IncLevel_d0r0",
                            "IncLevel2":"IncLevel_d1000r1","IncLevelDifference":"IncLevelDifference_d0r0_d1000r1"}, axis='columns')
#print(df_MBNL, df_FOX, df_MF)



#Sort and output total common events between rMATS datasets
#comapre d1000r0 to d0r1
df_MBNL_FOX = pandas.merge(df_MBNL, df_FOX, how = 'inner', on = ['GeneID', 'geneSymbol', 'strand', 'chr', 'exonStart_0base',
                                                                'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'])
df_MBNL_FOX = df_MBNL_FOX.drop(['ID_x', 'ID.1_x','ID_y','ID.1_y'], axis = 1)
print('df_MBNL_FOX\n', df_MBNL_FOX, '\n')

#compare d1000r0 to d1000r1
df_MBNL_MF = pandas.merge(df_MBNL, df_MF, how = 'inner', on = ['GeneID', 'geneSymbol', 'strand', 'chr', 'exonStart_0base',
                                                                'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'])
df_MBNL_MF = df_MBNL_MF.drop(['ID_x', 'ID.1_x','ID_y','ID.1_y'], axis = 1)
print('df_MBNL_MF\n', df_MBNL_MF, '\n')

#commpare d0r1 to d1000r1
df_FOX_MF = pandas.merge(df_FOX, df_MF, how = 'inner', on = ['GeneID', 'geneSymbol', 'strand', 'chr', 'exonStart_0base',
                                                            'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'])
df_FOX_MF = df_FOX_MF.drop(['ID_x', 'ID.1_x','ID_y','ID.1_y'], axis = 1)
print('df_FOX_MF\n', df_FOX_MF, '\n')

#compare all 3
df_comp_all = pandas.merge(df_FOX_MF, df_MBNL, how = 'inner', on = ['GeneID', 'geneSymbol', 'strand', 'chr', 'exonStart_0base',
                                                            'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'])
df_comp_all = df_comp_all.drop(['IJC_d0r0_y', 'SJC_d0r0_y','ID','ID.1', 'IJC_d0r0_y', 'SJC_d0r0_y'], axis = 1)
print('df_comp_all\n', df_comp_all, '\n')

#saved as ddCell_All_common_events
df_comp_all.to_csv('All_common_events.txt', index = None, sep = '\t', float_format = '%9.6f')


#abs() is to account for exclusion values(i.e. negative IncLevelDifference)
#sort_a = df[(df.FDR < 0.05) & (df.IncLevelDifference.abs() >= 0.1)]
#sort_b = sort_a.sort_values('FDR', ascending = True)
#print(sort_b)

#sort_b.to_excel('out_excel.xlsx', index = None)

#sort_b.to_csv('out.txt', index = None, sep = '\t', float_format = '%9.6f')
