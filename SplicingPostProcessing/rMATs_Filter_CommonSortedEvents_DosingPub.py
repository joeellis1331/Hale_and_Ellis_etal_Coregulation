#!/usr/bin/python3

import sys
import numpy #pandas is a dependeny of numpy
import pandas

#NOT .read_csv that uses a standard '' separator. read_table uses '\t'
df_MEF = pandas.read_table(sys.argv[1]) #reads MEF Common file
df_HEK = pandas.read_table(sys.argv[2]) #reads HEK commmon file
print('df_MEF\n',df_MEF, '\n', 'df_HEK\n', df_HEK, '\n',)

#Sort for MBNL significant regulated events (PSI > 0.1 and FDR < 0.05 and sorts with lowest FDR on top)
#abs() is to account for exclusion values(i.e. negative IncLevelDifference)
df_MEF = df_MEF[(df_MEF.FDR_d0r0_d1000r0 < 0.05) & (df_MEF.IncLevelDifference_d0r0_d1000r0.abs() >= 0.1)]
#Filter for read counts
df_MEF = df_MEF[((df_MEF.IJC_d0r0_x + df_MEF.SJC_d0r0_x) >= 5) &
                ((df_MEF.IJC_d0r1 + df_MEF.SJC_d0r1) >= 5) &
                ((df_MEF.IJC_d1000r1 + df_MEF.SJC_d1000r1) >= 5) &
                ((df_MEF.IJC_d1000r0 + df_MEF.SJC_d1000r0) >= 5)]
df_MEF = df_MEF.sort_values('IncLevelDifference_d0r0_d1000r0', ascending = True)
print(df_MEF)

df_HEK = df_HEK[(df_HEK.FDR_d0r0_d1000r0 < 0.05) & (df_HEK.IncLevelDifference_d0r0_d1000r0.abs() >= 0.1)]
df_HEK = df_HEK[((df_HEK.IJC_d0r0_x + df_HEK.SJC_d0r0_x) >= 5) &
                ((df_HEK.IJC_d0r1 + df_HEK.SJC_d0r1) >= 5) &
                ((df_HEK.IJC_d1000r1 + df_HEK.SJC_d1000r1) >= 5) &
                ((df_HEK.IJC_d1000r0 + df_HEK.SJC_d1000r0) >= 5)]
df_HEK = df_HEK.sort_values('IncLevelDifference_d0r0_d1000r0', ascending = True)
print(df_HEK)

#prints to excel file
df_MEF.to_excel('MBNLsig_common_ddMEF_FiltRead.xlsx', index = None)
df_HEK.to_excel('MBNLsig_common_ddHEK_FiltRead.xlsx', index = None)


#sort_b.to_csv('out.txt', index = None, sep = '\t', float_format = '%9.6f')
