# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 12:10:56 2013

@author: Daniel Struck
"""
from collections import Counter

zaehlen = Counter()
total = 0
with open("2013-08-30_pol_lanl.csv") as f:
    next(f)
    for line in f:
        token = line.split("\t")
        zaehlen[token[1].split(".")[0]] += 1
        total += 1

for e in sorted(zaehlen):
    print e.replace("_","\\textunderscore ") + " & " + str(zaehlen[e]) + " \\\\"
    
print total
        