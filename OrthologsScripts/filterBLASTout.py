#!/usr/bin/python3

import sys
import pandas
import re

#read results
df = pandas.read_excel("humanMouse_tBlastxResultsCP.xlsx")
print(df)

#sort %identity values from greatest to least
df = df.sort_values(by = ['%ident'], ascending = False)
print(df)

#remove duplicates keeping the first instance. Because I have sorted the values from high to low it keeps the greatest %ident of the duplicates
df = df.drop_duplicates(subset = ['qseqid', 'sseqid'], keep = 'first')
print(df)

#remove any that do not pass the %ident cutoff of 70%
df = df[df['%ident'] >= 70.0]
print(df)

#remove rows in which genes are not matching, blast aligned all the sequences to each other
for row in df.itertuples():
    que = re.match(r'.+?_', row.qseqid)
    sub = re.match(r'.+?_', row.sseqid)
    if que.group().lower() == sub.group().lower():
        pass
    else:
        df = df.drop(labels = row.Index)
print(df)

df.to_csv('conservedEventsFinal.csv', index = None)
